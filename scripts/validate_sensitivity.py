"""Compare the two submission sensitivity sheets with this run's fitted results."""
from pathlib import Path
import json, math
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/sensitivity'
REF=ROOT/'reference/submitted_results.xlsx'
ATOL=1e-8; RTOL=1e-7

def check(code,rows):
    sheets=pd.ExcelFile(REF).sheet_names
    sheet=next(s for s in sheets if s.startswith(code+'_'))
    old=pd.read_excel(REF,sheet_name=sheet)
    new=pd.DataFrame(rows,columns=old.columns)
    assert old.shape==new.shape
    mismatches=[];maxerr=0.;numbers=0
    for i in range(len(old)):
        for column in old.columns:
            a,b=new.loc[i,column],old.loc[i,column]
            if pd.isna(a) and pd.isna(b):continue
            if isinstance(a,(int,float,np.number)) and isinstance(b,(int,float,np.number)):
                numbers+=1;maxerr=max(maxerr,abs(float(a)-float(b)))
                same=math.isclose(float(a),float(b),rel_tol=RTOL,abs_tol=ATOL)
            else:same=str(a)==str(b)
            if not same:mismatches.append({'row':i+2,'column':column,'new':str(a),'reference':str(b)})
    new.to_csv(RUN/(code+'_regenerated_submission.tsv'),sep='\t',index=False,encoding='utf-8-sig')
    return {'sheet':sheet,'rows':len(new),'numeric_cells':numbers,'passed':not mismatches,
        'max_absolute_error':maxerr,'rtol':RTOL,'atol':ATOL,'mismatches':mismatches}

def main():
    n1=[]
    specs=[
        ('diversity/contrast_results.tsv','spec','with_diversity','原分类：含多样性','原2748人，40列饮食扩展基础设计'),
        ('diversity/contrast_results.tsv','spec','without_diversity','原分类：移除多样性','同2748人，移除Shannon及log丰富度'),
        ('provenance/round_contrasts.tsv','scheme','same_subset_baseline','轮次子集：未调轮次','同1386人，唯一匹配发表轮次'),
        ('provenance/round_contrasts.tsv','scheme','sample_count_pooled_round_fixed_effects','轮次子集：调整合并轮次','同1386人；14个n<20轮次合并为其他，共19层'),
        ('bloom/bloom_contrast_results.tsv','scenario','R07_baseline','R07：候选过滤前','同2748人，R07重新注释，含饮食PC'),
        ('bloom/bloom_contrast_results.tsv','scenario','R07_filtered','R07：候选过滤后','同2748人，去除全部7个精确匹配候选ASV后重建'),
        ('bloom/qc5000_contrast_results.tsv','scenario','R07_baseline','过滤后深度合格子集：过滤前','同2343人，以过滤后总读数≥5000固定子集'),
        ('bloom/qc5000_contrast_results.tsv','scenario','R07_filtered','过滤后深度合格子集：过滤后','同2343人，候选过滤后重建；不与全2748人直接比较')]
    for rel,key,value,label,note in specs:
        table=pd.read_csv(RUN/rel,sep='\t');table=table[table[key].eq(value)]
        if 'threshold' in table:table=table[table.threshold.eq(1)]
        assert len(table)==3 and table.valid.all()
        for r in table.to_dict('records'):
            q=r.get('q_BH_selected3',r.get('q_BH_3'))
            n1.append([label,r['background'],1,int(r['n']),r['ratio_OR'],r['low'],r['high'],r['p'],q,
                '乳酪杆菌属/乳杆菌属','同一方案3项BH；不替代原2340项q',note,'08_本轮方法复核/'+rel])
    table=pd.read_csv(RUN/'within_person/within_person_cluster_results.csv')
    assert len(table)==18 and table.valid.all()
    n2=[];names={'all':'全部信息样本','drop_most_sampled':'排除采样最多者'}
    for r in table.to_dict('records'):
        n2.append([int(r['threshold']),names[r['spec']],r['background'],r['target'],int(r['people']),int(r['samples']),int(r['events']),
            r['OR'],r['low'],r['high'],r['p'],r['q_BH_6'],r['cluster_se'],r['reference'],
            '每一阈值/规则6项BH','同期目标属内关联，非两目标属直接差异','08_本轮方法复核/within_person/within_person_cluster_results.csv'])
    checks=[check('N01',n1),check('N02',n2)]
    # Failure is a reproducible diagnostic, not a result to silently remove.
    rounds=pd.read_csv(RUN/'provenance/round_contrasts.tsv',sep='\t')
    invalid=rounds[rounds.scheme.eq('all_round_fixed_effects')]
    assert len(invalid)==3 and not invalid.valid.any()
    assert invalid[['ratio_OR','low','high','p','q_BH_3']].isna().all().all()
    report={'all_passed':all(x['passed'] for x in checks),'checks':checks,
        'invalid_full_round_contrasts_retained_as_missing':len(invalid),
        'new_fits_and_regenerated_tables':'runs/sensitivity',
        'scope':'N01 eight fixed scenarios and N02 eighteen participant-clustered conditional-logit fits; no external-cohort effect validation claimed.'}
    (ROOT/'reports/sensitivity_validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    if not report['all_passed']:raise SystemExit(1)

if __name__=='__main__':main()
