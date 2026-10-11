"""Bounded, resumable experiment. Run with BLAS threads limited to one.

No model selection receives validation or test coefficients. The truth is used
to generate observations and, after selection, to evaluate reconstructions.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time
import numpy as np
import scipy
from .data import design, datasets, noise, profiles
from .direct import analytic_matrix, fd_matrix, mask_indices
from .filters import fit_candidates, select_candidates, predict
from .metrics import metrics_for_predictions
from .storage import atomic_npz, write_array_chunks, validate_condition, bounded_disk_usage

METHODS = ["tikhonov_l2", "tikhonov_h1", "tsvd", "wiener", "wiener_cap",
           "learned_signed", "learned_signed_cap", "learned_positive",
           "learned_bounded", "learned_bounded_cap", "lmmse"]
METRICS = ["error_abs", "visible", "null", "residual_nominal", "residual_independent"]


def atomic_json(path, value):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    os.replace(temp, path)


def file_hash(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def cpu_time():
    return sum(resource.getrusage(k).ru_utime + resource.getrusage(k).ru_stime
               for k in [resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN])


def check_limits(output, started, previous_cpu, limit=14400):
    used = previous_cpu + cpu_time() - started
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    size = bounded_disk_usage(output.parent)
    if used >= limit or rss >= 2 * 1024**3 or size >= 200_000_000:
        raise RuntimeError("CPU, resident-memory or output-storage limit reached")
    return dict(cpu_seconds=used, peak_rss_bytes=rss, stored_bytes=size)


def run(mode, output, budget_cpu=14400, prerequisite=None, tests_report=None):
    cfg = design(mode)
    if tests_report is None:
        raise ValueError("Independent test report required")
    test_evidence = json.loads(Path(tests_report).read_text())
    if not test_evidence.get("success"):
        raise ValueError("Independent tests have not passed")
    for filename, expected in test_evidence["source_sha256"].items():
        if file_hash(Path(__file__).parent/filename) != expected:
            raise ValueError("Source changed after independent tests")
    prior_cpu = max(5., float(test_evidence["cpu_seconds"]))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if mode != "pilot":
        if prerequisite is None:
            raise ValueError("Campaign requires an explicit successful pilot report")
        pilot = json.loads(Path(prerequisite).read_text())
        prior_cpu = pilot["cpu_seconds"]
        if pilot.get("verdict") != "PASS":
            raise ValueError("Pilot has not passed")
        if pilot["feasibility"][mode]["estimated_total_cpu_seconds_with_margin"] > budget_cpu:
            raise RuntimeError("Pilot predicts budget exceeded for selected design")
    sources = {p.name: file_hash(p) for p in sorted(Path(__file__).parent.glob("*.py"))}
    config_record = dict(design=cfg, source_sha256=sources,
                         numpy=np.__version__, scipy=scipy.__version__,
                         python=sys.version, seed_root=20261011,
                         axes=["training", "mask", "noise", "regime", "object", "repeat", "method"],
                         methods=METHODS, metrics=METRICS,
                         entry_status={"success": 0, "numerical_failure": 1,
                                       "metric_failure": 2, "family_failure": 4, "missing": 128},
                         pilot_ids=dict(train=[0, 23], validation=[24, 31],
                                        test=[32, 39], noiseless=[40, 55]))
    cp = output / "configuration.json"
    if cp.exists() and json.loads(cp.read_text()) != config_record:
        raise RuntimeError("Resume configuration or source code differs; preserve existing run")
    if not cp.exists():
        atomic_json(cp, config_record)
    state_path = output / "state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else dict(completed=[], cpu_seconds=prior_cpu)
    if state.get("running"):
        # A process killed before its checkpoint has unknown CPU consumption.
        # Debit its full reservation conservatively; never silently reset it.
        state["cpu_seconds"] += state.get("reservation_cpu_seconds", budget_cpu)
    previous_cpu = state["cpu_seconds"]
    started, wall_started = cpu_time(), time.monotonic()
    limit = min(budget_cpu, 1800 if mode == "pilot" else 14400)
    if previous_cpu >= limit:
        raise RuntimeError("No remaining measured/conservatively reserved CPU budget")
    state.update(running=True, reservation_cpu_seconds=min(600.,limit-previous_cpu) if mode=="pilot" else limit-previous_cpu)
    atomic_json(state_path, state)
    shape = (cfg["games"], len(cfg["masks"]), len(cfg["levels"]),
             len(cfg["regimes"]), cfg["test"], cfg["repetitions"], 11)
    expected_storage = int(np.prod(shape)) * 41 + 40_000_000
    if expected_storage > 200_000_000:
        raise RuntimeError("Preflight output allocation exceeds fixed budget")
    existing_bytes = sum((output/(name+".npy")).stat().st_size for name in METRICS+["status"] if (output/(name+".npy")).exists())
    bounded_disk_usage(output.parent, planned_bytes=max(0,expected_storage-existing_bytes))
    arrays = {}
    for name in METRICS + ["status"]:
        path = output / (name + ".npy")
        exists = path.exists()
        arrays[name] = np.lib.format.open_memmap(path, mode="r+" if exists else "w+",
                     dtype=np.uint8 if name == "status" else np.float64, shape=shape)
        if arrays[name].shape != shape or arrays[name].dtype != np.dtype(np.uint8 if name=="status" else np.float64):
            raise RuntimeError("Corrupt resume array shape or dtype")
        if not exists:
            arrays[name][:] = 128 if name == "status" else np.nan
            arrays[name].flush()
    train, val, nominal, rough, vi, ti = datasets(cfg)
    for key in state["completed"]:
        validate_condition(arrays, tuple(map(int,key.split("_"))), state["condition_hashes"][key])
    truths = {f"train_{g}": c for g,c in enumerate(train)}
    truths.update(validation=val, nominal=nominal, rough=rough if rough is not None else np.empty((0,cfg["K"])))
    truth_manifest={}
    for name,values in truths.items():
        chunks=write_array_chunks(output/"truth"/name,values)
        for chunk in chunks: chunk["path"]=str(Path(chunk["path"]).relative_to(output))
        truth_manifest[name]=chunks
    atomic_json(output/"truth_manifest.json",truth_manifest)
    control_report = {}
    if pilot_mode := (mode == "pilot"):
        tfine=cpu_time()
        fine = fd_matrix(9, 256, 64)
        exact = analytic_matrix(9)
        relative = np.linalg.norm(fine-exact,axis=0)/np.linalg.norm(exact,axis=0)
        if relative[0] > 1e-10 or np.max(relative[1:]) > .005:
            raise RuntimeError("Fine-grid modal control failed")
        controls = profiles(0, 0, range(40,56), cfg["K"])
        exact_controls = controls @ analytic_matrix(cfg["K"]).T
        if not np.isfinite(exact_controls).all():
            raise RuntimeError("Nonfinite noiseless control")
        control_report=dict(fine_grid_cpu_seconds=cpu_time()-tfine,
                            modal_relative_errors=relative.tolist(),
                            noiseless_object_ids=list(range(40,56)))
        atomic_npz(output/"noiseless_controls.npz",coefficients=controls,observations=exact_controls)
    inverse = fd_matrix(cfg["n"], cfg["nx"], cfg["ny"])
    generators = {"nominal": analytic_matrix(cfg["K"]),
                  "rough": analytic_matrix(cfg["K"]),
                  "height_2": analytic_matrix(cfg["K"], H=.25*1.02),
                  "height_5": analytic_matrix(cfg["K"], H=.25*1.05)}
    pilot_mode = mode == "pilot"
    val_z = noise(0 if pilot_mode else 2, 0, vi, [0, 1])
    test_z = {"nominal": noise(0 if pilot_mode else 3, 0, ti, range(cfg["repetitions"]+1))}
    if not pilot_mode:
        test_z["rough"] = noise(4, 0, ti, range(cfg["repetitions"]+1))
    timing = state.get("timing", dict(fit_cpu=0., test_cpu=0., conditions=0, reconstructions=0))
    failures = 0
    samples = {}
    for sample_path in sorted(output.glob("failure_samples_*.npz")):
        with np.load(sample_path) as saved:
            samples.update({tuple(map(int,key)): value for key,value in zip(saved["ids"],saved["values"])})
    for g in range(cfg["games"]):
        clean_train = train[g] @ generators["nominal"].T
        scale = float(np.sqrt(np.mean(clean_train**2)))
        d0 = float(1e-12 * np.median(np.linalg.norm(train[g], axis=1)))
        train_z = noise(0 if pilot_mode else 1, g, range(cfg["train"]), [0])[:, 0]
        for mi, mask in enumerate(cfg["masks"]):
            ix = mask_indices(mask)
            weighted = inverse[ix] / 8
            a = float(np.linalg.norm(weighted, 2))
            qscale = 1/(8*a)
            A = inverse[ix] * qscale
            _, singular, vt = np.linalg.svd(A, full_matrices=False)
            rank = int(np.count_nonzero(singular > 1e-12*singular[0]))
            Vr = vt[:rank].T
            for ni, nu in enumerate(cfg["levels"]):
                key = f"{g}_{mi}_{ni}"
                if key in state["completed"]:
                    continue
                check_limits(output, started, previous_cpu, limit)
                if pilot_mode and time.monotonic()-wall_started >= 600:
                    raise RuntimeError("Pilot wall-time limit reached")
                tau = scale*nu
                dtrain = (clean_train + tau*train_z)[:, ix] * qscale
                vclean = val @ generators["nominal"].T
                dval = (vclean[:, None, :] + tau*val_z)[:, :, ix] * qscale
                tfit = cpu_time()
                candidates = fit_candidates(A, train[g][:, :cfg["n"]], dtrain,
                                             (tau*qscale)**2)
                selected, diagnostics = select_candidates(candidates, A, dval[:, 0], dval[:, 1], return_diagnostics=True)
                for diagnostic in diagnostics:
                    for field,value in diagnostic.items():
                        if isinstance(value,float) and not np.isfinite(value): diagnostic[field]=None
                atomic_json(output/("candidates_"+key+".json"),diagnostics)
                fit_elapsed = cpu_time()-tfit
                if set(selected) != set(METHODS):
                    raise RuntimeError("Incomplete method selection")
                invalid = [m for m in METHODS if not selected[m]["valid"]]
                if pilot_mode and invalid:
                    raise RuntimeError("Invalid pilot families: "+str(invalid))
                payload = {"A": A, "singular": singular, "Vr": Vr,
                           "scale": np.array(scale), "d0": np.array(d0), "qscale": np.array(qscale)}
                summary = {}
                for m in METHODS:
                    e = selected[m]
                    for field in ["R", "mu_c", "mu_d"]:
                        payload[m+"__"+field] = e[field]
                    summary[m] = {k: (float(v) if np.isscalar(v) and not isinstance(v, (str, bool)) else v)
                                  for k,v in e.items() if k in ["parameter", "grid_index", "operator_norm", "valid", "score"]}
                    for k,v in summary[m].items():
                        if isinstance(v, float) and not np.isfinite(v): summary[m][k] = None
                atomic_npz(output / ("estimators_"+key+".npz"), **payload)
                atomic_json(output / ("selection_"+key+".json"), summary)
                ttest = cpu_time()
                for ri, regime in enumerate(cfg["regimes"]):
                    truth = rough if regime == "rough" else nominal
                    clean = truth @ generators[regime].T
                    z = test_z["rough" if regime == "rough" else "nominal"]
                    d = (clean[:, None, :] + tau*z[:, :cfg["repetitions"]])[:, :, ix]*qscale
                    d2 = (clean + tau*z[:, cfg["repetitions"]])[:, ix]*qscale
                    flat_d = d.reshape(-1, len(ix))
                    flat_truth = np.repeat(truth, cfg["repetitions"], axis=0)
                    independent = np.repeat(d2, cfg["repetitions"], axis=0)
                    for method_index, m in enumerate(METHODS):
                        chat = predict(selected[m], flat_d)
                        result = metrics_for_predictions(chat, flat_truth, A, flat_d, independent, Vr)
                        if not selected[m]["valid"]:
                            result["status"][:] |= 4
                            for metric in METRICS: result[metric][:]=np.inf
                        denominator=np.maximum(np.linalg.norm(flat_truth,axis=1),d0)
                        with np.errstate(over="ignore",invalid="ignore"):
                            bad_relative=~np.isfinite(result["error_abs"]/denominator)
                        newly_failed=bad_relative & (result["status"]==0)
                        result["status"][newly_failed] |= 2
                        for metric in METRICS:
                            result[metric][newly_failed]=np.inf
                        condition_id=((g*len(cfg["masks"])+mi)*len(cfg["levels"])+ni)*len(cfg["regimes"])+ri
                        for row in np.flatnonzero(result["status"]):
                            sid=(1,condition_id,method_index,int(row))
                            if len(samples)<1024 or sid<max(samples):
                                raw=np.full(128,np.nan)
                                raw[:cfg["n"]]=chat[row]
                                raw[64:64+len(ix)]=flat_d[row]
                                samples[sid]=raw
                                if len(samples)>1024: del samples[max(samples)]
                        for name in METRICS+["status"]:
                            arrays[name][g,mi,ni,ri,:,:,method_index] = result[name].reshape(cfg["test"],cfg["repetitions"])
                        failures += int(np.count_nonzero(result["status"]))
                for arr in arrays.values(): arr.flush()
                count = len(cfg["regimes"])*cfg["test"]*cfg["repetitions"]*11
                timing["fit_cpu"] += fit_elapsed
                timing["test_cpu"] += cpu_time()-ttest
                timing["conditions"] += 1
                timing["reconstructions"] += count
                state.setdefault("condition_hashes",{})[key]=validate_condition(arrays,(g,mi,ni))
                state["completed"].append(key)
                state.update(check_limits(output, started, previous_cpu, limit))
                state["timing"] = timing
                if samples:
                    ids=sorted(samples)
                    for chunk,start in enumerate(range(0,len(ids),512)):
                        group=ids[start:start+512]
                        atomic_npz(output/f"failure_samples_{chunk}.npz",ids=np.asarray(group,dtype=np.uint64),values=np.asarray([samples[i] for i in group]))
                if pilot_mode and time.monotonic()-wall_started > 600:
                    atomic_json(state_path,state)
                    raise RuntimeError("Pilot exceeded600seconds; no PASS")
                atomic_json(state_path, state)
                print(json.dumps(dict(condition=key, cpu=state["cpu_seconds"], failed=failures)),flush=True)
    complete = not np.any(arrays["status"] == 128)
    total_failures=int(np.count_nonzero(arrays["status"] & 7))
    report = dict(mode=mode, complete=bool(complete), verdict="PASS" if complete and (not pilot_mode or total_failures==0) else "FAIL",
                  failures=int(np.count_nonzero(arrays["status"] & 7)),
                  **check_limits(output, started, previous_cpu, limit),
                  elapsed_seconds=time.monotonic()-wall_started, timing=timing,
                  reconstruction_count=int(np.prod(shape)), controls=control_report,
                  independent_tests_report_sha256=file_hash(tests_report),
                  prior_cpu_debit_seconds=prior_cpu)
    if pilot_mode:
        feasibility={}
        for target in ["full","fallback"]:
            c=design(target)
            # Cubic dimension scaling conservatively covers dense factorizations;
            # sample scaling covers moment accumulation and predictive validation.
            dims=(c["n"]/cfg["n"])**3
            fit=timing["fit_cpu"]/timing["conditions"]*27*max(c["train"]/24,c["validation"]/8)*dims
            count=27*len(c["regimes"])*c["test"]*c["repetitions"]*11
            test=timing["test_cpu"]/timing["reconstructions"]*count*(c["n"]/cfg["n"])**2
            estimate=report["cpu_seconds"]+1.5*(fit+test+1800+1800)
            feasibility[target]=dict(estimated_total_cpu_seconds_with_margin=estimate,
                                    fit_cpu=fit,test_cpu=test,preparation_and_analysis_reserve=3600)
        report["feasibility"]=feasibility
    state.update(running=False,reservation_cpu_seconds=0,cpu_seconds=report["cpu_seconds"])
    atomic_json(state_path,state)
    atomic_json(output/"report.json",report)
    manifest={str(p.relative_to(output)):dict(bytes=p.stat().st_size,sha256=file_hash(p))
              for p in sorted(output.rglob("*")) if p.is_file() and p.name!="manifest.json"}
    atomic_json(output/"manifest.json",manifest)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode",choices=["pilot","full","fallback"],required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--pilot-report")
    p.add_argument("--tests-report",required=True)
    args=p.parse_args()
    print(json.dumps(run(args.mode,args.output,prerequisite=args.pilot_report,tests_report=args.tests_report),indent=2))


if __name__ == "__main__":
    main()
