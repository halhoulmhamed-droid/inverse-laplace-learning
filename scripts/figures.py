"""Frozen-rule figures from completed outputs; never fits or draws new seeds.

Run only after analysis, with its JSON report and cumulative CPU debit.
PDF/PNG/NPZ plus manifest are bounded jointly by 10,000,000 bytes.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from laplace_learning.data import noise
from laplace_learning.direct import analytic_matrix, mask_indices
from laplace_learning.metrics import ordered_quantile

LABELS = ['Tikhonov L2', 'Tikhonov H1', 'TSVD', 'Wiener', 'Wiener plafonné',
          'Appris signé', 'Appris signé plafonné', 'Appris positif',
          'Appris borné', 'Appris borné plafonné', 'LMMSE']
MASKS = ['Complet', 'Alterné', 'Contigu']
BASE_URL = 'https://github.com/halhoulmhamed-droid/inverse-laplace-learning'


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def cpu():
    return sum(resource.getrusage(k).ru_utime + resource.getrusage(k).ru_stime
               for k in [resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN])


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def link_manifest(directory, commit):
    """Bind existing artifacts only to a verified local Git commit containing them."""
    if len(commit) != 40 or any(c not in '0123456789abcdef' for c in commit):
        raise ValueError('Supply an exact lowercase 40-character commit SHA')
    path = directory/'manifest.json'
    manifest = json.loads(path.read_text())
    for entry in manifest['figures']:
        links = {}
        for kind, rel in entry['paths'].items():
            actual = subprocess.check_output(['git', 'show', f'{commit}:{rel}'], cwd=ROOT)
            if hashlib.sha256(actual).hexdigest() != entry['sha256'][kind]:
                raise ValueError('Git artifact differs: '+rel)
            links[kind] = f'{BASE_URL}/blob/{commit}/{rel}'
        entry['github_links'] = links
    manifest['git_commit'] = commit
    write_json(path, manifest)


def run(args):
    started = 0.0  # Account for imports and all process CPU, not rendering alone.
    source, out = Path(args.results).resolve(), Path(args.output).resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError('Preserve previous figures: output must be absent or empty')
    config_path = source/'configuration.json'
    config = json.loads(config_path.read_text())
    report = json.loads((source/'report.json').read_text())
    analysis = json.loads(Path(args.analysis_report).read_text())
    if not report.get('complete') or report.get('verdict') != 'PASS' or not analysis.get('complete'):
        raise ValueError('Completed campaign and analysis required')
    if args.prior_cpu_seconds < report['cpu_seconds']:
        raise ValueError('Cumulative CPU debit cannot precede the campaign')
    for module in ['data.py', 'direct.py', 'metrics.py']:
        if sha(ROOT/'laplace_learning'/module) != config['source_sha256'][module]:
            raise ValueError('Reconstruction module changed since campaign: '+module)
    cfg = config['design']
    if cfg['mode'] != 'full' or len(config['methods']) != 11:
        raise ValueError('This figure specification requires the full fixed design')
    source_manifest = json.loads((source/'manifest.json').read_text())
    used_sources = {}
    def checked(rel):
        path = source/rel
        observed = sha(path)
        if source_manifest[rel]['sha256'] != observed:
            raise ValueError('Source hash mismatch: '+rel)
        used_sources[str(path.relative_to(ROOT))] = observed
        return path
    for rel in ['configuration.json', 'report.json', 'truth_manifest.json']:
        checked(rel)
    metrics = {name: np.load(checked(name+'.npy'), mmap_mode='r')
               for name in ['error_abs', 'residual_nominal', 'status']}
    if np.any(metrics['status'] == 128):
        raise ValueError('Missing experiment rows')
    tm = json.loads((source/'truth_manifest.json').read_text())
    truth = np.concatenate([np.load(checked(c['path'])) for c in tm['nominal']])
    norms = np.linalg.norm(truth, axis=1)
    estimators = {}
    def estimator(g, mi):
        key = (g, mi)
        if key not in estimators:
            with np.load(checked(f'estimators_{g}_{mi}_1.npz')) as z:
                estimators[key] = dict(z)
        return estimators[key]
    def errors(g, mi, ri):
        return metrics['error_abs'][g, mi, 1, ri] / np.maximum(norms, estimator(g,mi)['d0'])[:,None,None]
    out.mkdir(parents=True)
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9,
                         'axes.spines.top':False, 'axes.spines.right':False,
                         'pdf.fonttype':42, 'savefig.dpi':120})
    manifest = dict(figures=[], git_commit=None, trace_selections=[], absences=[],
                    source_sha256=used_sources, analysis_report_sha256=sha(args.analysis_report),
                    scope='All training games shown separately; repeated objects are not independent training replicates.',
                    figure_budget_bytes=10_000_000)
    def check():
        size = sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
        elapsed = cpu()-started
        if size >= 10_000_000 or elapsed > 1800 or elapsed+args.prior_cpu_seconds > 14400:
            raise RuntimeError('Figure storage or analysis CPU limit exceeded; preserve partial outputs')
    def save(name, fig, data, description):
        check()
        paths = {'data':out/(name+'.npz'), 'pdf':out/(name+'.pdf'), 'png':out/(name+'.png'),
                 'script':Path(__file__).resolve(), 'configuration':config_path}
        np.savez_compressed(paths['data'], **data)
        fig.savefig(paths['pdf'], bbox_inches='tight', metadata={'Title':description})
        fig.savefig(paths['png'], bbox_inches='tight')
        plt.close(fig)
        manifest['figures'].append(dict(name=name, description=description,
              paths={k:str(v.relative_to(ROOT)) for k,v in paths.items()},
              sha256={k:sha(v) for k,v in paths.items()}, github_links=None))
        check()
    colors = list(plt.get_cmap('tab20').colors[:11])
    fig, ax = plt.subplots(figsize=(7,4))
    spectra = {}
    for mi in range(3):
        s = estimator(0,mi)['singular']; spectra[f'mask_{mi}'] = s
        ax.semilogy(np.arange(1,len(s)+1),s,label=MASKS[mi])
    ax.axhline(1e-12,color='gray',ls=':',label='Seuil relatif 10⁻¹²')
    ax.set(xlabel='Indice singulier',ylabel='Valeur singulière normalisée');ax.legend()
    save('spectres',fig,spectra,'Spectres des opérateurs inverses, trois masques')
    summary = np.empty((3,3,2,11))
    for g in range(3):
        for mi in range(3):
            for qi,q in enumerate([.5,.95]):
                for j in range(11): summary[g,mi,qi,j]=ordered_quantile(errors(g,mi,0)[:,:,j],q)
    fig, axes = plt.subplots(2,3,figsize=(12,8),sharex=True)
    for mi in range(3):
        for qi in range(2):
            ax=axes[qi,mi]
            for g in range(3): ax.plot(range(11),summary[g,mi,qi],marker='o',ms=3,label=f'Entraînement {g}')
            ax.set_yscale('log');ax.set_title(MASKS[mi]+(' — médiane' if qi==0 else ' — q95'))
            ax.set_xticks(range(11),LABELS,rotation=85);ax.set_ylabel('Erreur relative L2')
    axes[0,0].legend();fig.tight_layout()
    save('comparaison_nominale',fig,{'quantiles':summary,'methods':np.array(LABELS)},'Toutes les méthodes, bruit 1 %, nominal ; trois entraînements séparés')
    heights = np.array([0,2,5]);values=np.empty((3,3,11))
    fig,axes=plt.subplots(1,3,figsize=(13,4),sharey=True)
    for g,ax in enumerate(axes):
        for j in range(11):
            for h,ri in enumerate([0,2,3]):values[g,h,j]=ordered_quantile(errors(g,0,ri)[:,:,j],.95)
            ax.plot(heights,values[g,:,j],marker='o',color=colors[j],ls=['-','--',':'][j%3],label=LABELS[j])
        ax.set(xlabel='Erreur de hauteur (%)',title=f'Entraînement {g}',yscale='log',xticks=heights)
    axes[0].set_ylabel('q95 de l’erreur relative L2');fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=4,bbox_to_anchor=(.5,-.16))
    save('effet_hauteur',fig,{'height_percent':heights,'q95':values},'Masque complet, bruit 1 % ; hauteur vraie modifiée sans réentraînement')
    contrasts=[('P1',2,0),('P2',0,3)]
    ratios=np.empty((2,3,2));fig,axes=plt.subplots(1,2,figsize=(9,4))
    for ci,(name,mi,ri) in enumerate(contrasts):
        for g in range(3):
            for qi,q in enumerate([.5,.95]):
                source_error=errors(g,mi,0 if qi==0 else ri)
                den=ordered_quantile(source_error[:,:,5],q);num=ordered_quantile(source_error[:,:,6],q)
                ratios[ci,g,qi]=num/den if np.isfinite(den) and den>0 else np.nan
        for qi,label in enumerate(['Médiane nominale','q95 contraste']):axes[ci].plot(range(3),ratios[ci,:,qi],marker='o',label=label)
        axes[ci].axhline(.90,color='gray',ls=':');axes[ci].axhline(1.05,color='gray',ls='--')
        axes[ci].set(title=name,xlabel='Entraînement',ylabel='Ratio plafonné / signé',xticks=[0,1,2]);axes[ci].legend()
    save('dispersion_entrainements',fig,{'ratios':ratios},'Ratios ponctuels par entraînement ; aucune couverture simultanée revendiquée')
    fig,axes=plt.subplots(1,2,figsize=(10,4));scatter_data={}
    for ci,(name,mi,ri) in enumerate(contrasts):
        e=errors(0,mi,ri).mean(axis=1);r=metrics['residual_nominal'][0,mi,1,ri].mean(axis=1)
        scatter_data[name+'_errors']=e;scatter_data[name+'_residuals']=r
        for j in range(11):axes[ci].scatter(r[:,j],e[:,j],s=5,alpha=.45,color=colors[j],label=LABELS[j],rasterized=True)
        axes[ci].set(xscale='linear',yscale='log',title=name+' — entraînement 0',xlabel='Résidu nominal moyen',ylabel='Erreur relative moyenne')
        axes[ci].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
        hidden=int(np.count_nonzero(~np.isfinite(r)|~np.isfinite(e)|(r<=0)|(e<=0)))
        axes[ci].text(.02,.02,f'{hidden} couples non représentables (conservés)',transform=axes[ci].transAxes,fontsize=6)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=4,bbox_to_anchor=(.5,-.17))
    save('residu_erreur',fig,scatter_data,'Moyennes objet sur dix réalisations ; toutes méthodes, aucun filtrage succès seuls')
    x=np.linspace(0,1,257);basis=np.sqrt(2)*np.cos(np.pi*x[:,None]*np.arange(cfg['K']));basis[:,0]=1
    def trace(name,mi,ri,obj,rep,reason):
        es=estimator(0,mi);height=.25*(1.05 if ri==3 else 1)
        z=noise(3,0,[obj],[rep])[0,0]
        d=(analytic_matrix(cfg['K'],H=height)@truth[obj]+cfg['levels'][1]*float(es['scale'])*z)[mask_indices(cfg['masks'][mi])]*float(es['qscale'])
        predictions=[];fig,axes=plt.subplots(3,4,figsize=(12,8),sharex=True)
        actual= basis@truth[obj];status=metrics['status'][0,mi,1,ri,obj,rep]
        for j,ax in enumerate(axes.flat):
            ax.plot(x,actual,color='black',lw=1,label='Vérité')
            if j==11:
                ax.clear();ax.axis('off')
                ax.text(.05,.8,'Noir : vérité\nCouleur : reconstruction\nÉchelles verticales propres\nà chaque méthode',transform=ax.transAxes,va='top')
                continue
            method=config['methods'][j];c=es[method+'__mu_c']+es[method+'__R']@(d-es[method+'__mu_d']);predictions.append(c)
            if status[j]==0 and np.isfinite(c).all():ax.plot(x,basis[:,:cfg['n']]@c,color=colors[j],lw=1,label='Reconstruction')
            else:ax.text(.05,.5,f'Sortie en échec (statut {int(status[j])})',transform=ax.transAxes)
            ax.set_title(LABELS[j]);ax.set_xlabel('x');ax.set_ylabel('q(x)')
        emean=errors(0,mi,ri)[obj,:,5].mean();rmean=metrics['residual_nominal'][0,mi,1,ri,obj,:,5].mean()
        fig.suptitle(f'{name} — objet {obj}, réalisation {rep} ; {reason}\nSigné : erreur moyenne={emean:.4g}, résidu moyen={rmean:.4g}')
        fig.tight_layout()
        manifest['trace_selections'].append(dict(name=name,object_id=int(obj),repetition=int(rep),reason=reason))
        save(name,fig,dict(x=x,truth_coefficients=truth[obj],reconstructed_coefficients=np.array(predictions),observations=d,status=status),reason)
    for contrast,mi,ri in contrasts:
        means=errors(0,mi,ri)[:,:,5].mean(axis=1);res=metrics['residual_nominal'][0,mi,1,ri,:,:,5].mean(axis=1)
        ordered=np.lexsort((np.arange(len(means)),means))
        for label,q in [('mediane',.5),('q95',.95),('pire',1.)]:
            obj=int(ordered[math.ceil(q*len(means))-1]);trace(contrast+'_'+label,mi,ri,obj,0,label+' selon erreur moyenne du filtre signé')
        low=ordered_quantile(res,.1);high=ordered_quantile(means,.9)
        ids=np.flatnonzero(np.isfinite(res)&np.isfinite(means)&(res<=low)&(means>=high))
        if len(ids):
            obj=int(sorted(ids,key=lambda i:(-means[i],res[i],i))[0]);trace(contrast+'_faible_residu',mi,ri,obj,0,'Déciles préfixés : faible résidu / erreur élevée')
        else:manifest['absences'].append(contrast+': aucun objet éligible faible résidu / erreur élevée')
        failed=np.argwhere(np.any(metrics['status'][0,mi,1,ri]!=0,axis=-1))
        if len(failed):trace(contrast+'_echec',mi,ri,int(failed[0,0]),int(failed[0,1]),'Premier échec lexicographique')
        else:manifest['absences'].append(contrast+': aucun échec numérique')
    manifest.update(cpu_seconds_this_run=cpu()-started, prior_cpu_seconds=args.prior_cpu_seconds,
                    cpu_seconds_cumulative=args.prior_cpu_seconds+cpu()-started,
                    complete=True, matplotlib=matplotlib.__version__, numpy=np.__version__)
    write_json(out/'manifest.json',manifest);check()
    print(json.dumps({'figures':len(manifest['figures']),'cpu_seconds_cumulative':manifest['cpu_seconds_cumulative'],'absences':manifest['absences']},ensure_ascii=False))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',default=str(ROOT/'results/full'))
    p.add_argument('--output',default=str(ROOT/'figures'))
    p.add_argument('--analysis-report')
    p.add_argument('--prior-cpu-seconds',type=float)
    p.add_argument('--link-commit')
    a=p.parse_args()
    if a.link_commit:link_manifest(Path(a.output),a.link_commit)
    else:
        if a.analysis_report is None or a.prior_cpu_seconds is None:p.error('Analysis report and cumulative CPU debit required')
        run(a)
