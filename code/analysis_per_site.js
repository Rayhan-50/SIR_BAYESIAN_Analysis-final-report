#!/usr/bin/env node
// ============================================================
// Bayesian SIR Model — Per Sentinel Site Analysis (Node.js)
// Raihan Research Group | 2026
//
// MCMC: Metropolis-within-Gibbs (per-parameter updates)
//       4 chains × 6000 iter (3000 warmup)
//       Per-parameter adaptive scaling → target 44% acceptance
// Plots: Interactive Plotly HTML
// ============================================================
'use strict';

const fs   = require('fs');
const path = require('path');

const PROJ = path.resolve(__dirname, '..');
const OUT  = path.join(PROJ, 'outputs');
fs.mkdirSync(OUT, { recursive: true });

// ── Site definitions ──────────────────────────────────────────
// SARI: hospital surveillance — use N=40000 (ward catchment)
//       and I0=50 (pre-existing infected at surveillance start)
const SITES = [
  { label:'ILI',    csv:'site_ILI.csv',    pop:11750,  I0:1,  desc:'ILI Outpatient Sentinel' },
  { label:'Severe', csv:'site_Severe.csv',  pop:11750,  I0:1,  desc:'Severe/Hospital Sentinel' },
  { label:'SARI',   csv:'site_SARI.csv',    pop:40000,  I0:50, desc:'SARI Sentinel' }
];

const N_CHAINS = 4;
const N_ITER   = 6000;
const N_WARMUP = 3000;

// Parameter names in log-space order
const PNAMES = ['logR0','logGamma','logRho','logPhi'];

// ═══════════════════════════════════════════════════════════════
// UTILITIES
// ═══════════════════════════════════════════════════════════════
function makeRng(seed) {
  let s = (seed >>> 0) || 1;
  return {
    rand()  { s^=s<<13; s^=s>>>17; s^=s<<5; return (s>>>0)/4294967296; },
    randn() {
      const u=Math.max(this.rand(),1e-15);
      return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*this.rand());
    }
  };
}

function readCSV(file) {
  const rows = fs.readFileSync(file,'utf8').split(/\r?\n/).filter(Boolean);
  const hdr  = rows[0].split(',').map(h=>h.trim());
  return rows.slice(1)
    .map(r=>{ const v=r.split(','); const o={}; hdr.forEach((h,i)=>o[h]=v[i]?v[i].trim():''); return o; })
    .filter(r=>r.cases!=='');
}

// ═══════════════════════════════════════════════════════════════
// SIR ODE (RK4) — accepts custom I0
// ═══════════════════════════════════════════════════════════════
function solveSIR(beta, gamma, N, T, I0=1) {
  let S=N-I0, I=I0, R=0;
  const It=[];
  for (let t=0; t<T; t++) {
    const dS=-beta*S*I/N, dI=beta*S*I/N-gamma*I;
    const k1={ds:dS,di:dI};
    const k2={ds:-beta*(S+.5*k1.ds)*(I+.5*k1.di)/N,
              di: beta*(S+.5*k1.ds)*(I+.5*k1.di)/N-gamma*(I+.5*k1.di)};
    const k3={ds:-beta*(S+.5*k2.ds)*(I+.5*k2.di)/N,
              di: beta*(S+.5*k2.ds)*(I+.5*k2.di)/N-gamma*(I+.5*k2.di)};
    const k4={ds:-beta*(S+k3.ds)*(I+k3.di)/N,
              di: beta*(S+k3.ds)*(I+k3.di)/N-gamma*(I+k3.di)};
    S+=(k1.ds+2*k2.ds+2*k3.ds+k4.ds)/6;
    I+=(k1.di+2*k2.di+2*k3.di+k4.di)/6;
    S=Math.max(0,S); I=Math.max(1e-9,I);
    It.push(I);
  }
  return It;
}

// ═══════════════════════════════════════════════════════════════
// LOG-GAMMA (Lanczos)
// ═══════════════════════════════════════════════════════════════
const LC=[0.99999999999980993,676.5203681218851,-1259.1392167224028,
          771.32342877765313,-176.61502916214059,12.507343278686905,
          -0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
function lgamma(z){
  if(z<0.5) return Math.log(Math.PI)-Math.log(Math.sin(Math.PI*z))-lgamma(1-z);
  z--; let x=LC[0]; for(let i=1;i<9;i++) x+=LC[i]/(z+i);
  const t=z+7.5;
  return .5*Math.log(2*Math.PI)+(z+.5)*Math.log(t)-t+Math.log(x);
}

function nbLogPMF(k,mu,phi){
  mu=Math.max(mu,1e-10);
  if(!isFinite(mu)||phi<=0) return -Infinity;
  return lgamma(phi+k)-lgamma(phi)-lgamma(k+1)
        +phi*Math.log(phi/(phi+mu))+k*Math.log(mu/(phi+mu));
}

// ═══════════════════════════════════════════════════════════════
// LOG-POSTERIOR
// theta = [logR0, logGamma, logRho, logPhi]
// ═══════════════════════════════════════════════════════════════
function logPost(theta, cases, N, I0) {
  const [logR0,logG,logRho,logPhi]=theta;
  const R0=Math.exp(logR0), gamma=Math.exp(logG),
        rho=Math.exp(logRho), phi=Math.exp(logPhi);
  const beta=R0*gamma;

  if(R0<=0||gamma<=0||rho<=0||rho>=1||phi<=0) return -Infinity;
  if(beta>30||gamma>10) return -Infinity;

  // Log-priors (with Jacobian for log-transform)
  const lp = (x,mu,sig) => -0.5*((x-mu)/sig)**2-Math.log(sig)-0.5*Math.log(2*Math.PI);
  const lpR0  = lp(logR0, Math.log(1.7), 0.25) - logR0;
  const lpG   = lp(logG,  Math.log(0.44),0.30) - logG;
  const lpRho = lp(logRho,Math.log(0.05),0.50) - logRho;
  const lpPhi = -phi - logPhi;  // Exponential(1) + Jacobian

  let It;
  try { It=solveSIR(beta,gamma,N,cases.length,I0); }
  catch{ return -Infinity; }

  let ll=0;
  for(let t=0;t<cases.length;t++){
    const mu=Math.max(rho*It[t],1e-10);
    const lp_=nbLogPMF(cases[t],mu,phi);
    if(!isFinite(lp_)) return -Infinity;
    ll+=lp_;
  }
  return lpR0+lpG+lpRho+lpPhi+ll;
}

// ═══════════════════════════════════════════════════════════════
// METROPOLIS-WITHIN-GIBBS (one parameter at a time)
// Each parameter gets its own adaptive scale → target 44%
// ═══════════════════════════════════════════════════════════════
function runChain(cases, N, I0, seed) {
  const rng=makeRng(seed);

  // Starting values: perturb gently around prior mode
  let cur=[
    Math.log(1.7) +0.05*rng.randn(),
    Math.log(0.44)+0.05*rng.randn(),
    Math.log(0.05)+0.08*rng.randn(),
    Math.log(1.0) +0.10*rng.randn()
  ];
  let curLP=logPost(cur,cases,N,I0);
  // If start is invalid, nudge to prior mode
  if(!isFinite(curLP)){
    cur=[Math.log(1.7),Math.log(0.44),Math.log(0.05),Math.log(1.0)];
    curLP=logPost(cur,cases,N,I0);
  }

  // Per-parameter proposal scales and acceptance counts
  let scales=[0.05,0.05,0.08,0.12];
  let acc=[0,0,0,0], att=[0,0,0,0];

  const draws=[];

  for(let iter=0;iter<N_ITER;iter++){
    // Update each parameter independently (Metropolis-within-Gibbs)
    for(let pi=0;pi<4;pi++){
      const prop=[...cur];
      prop[pi]+=scales[pi]*rng.randn();
      const propLP=logPost(prop,cases,N,I0);
      att[pi]++;
      if(Math.log(Math.max(rng.rand(),1e-300))<propLP-curLP){
        cur=prop; curLP=propLP; acc[pi]++;
      }
    }

    // Adapt per-parameter scales every 50 steps during warmup
    if(iter<N_WARMUP&&(iter+1)%50===0){
      for(let pi=0;pi<4;pi++){
        const rate=acc[pi]/att[pi];
        const f=rate>0.50?1.20:rate>0.44?1.05:rate<0.20?0.80:rate<0.30?0.95:1.0;
        scales[pi]=Math.min(Math.max(scales[pi]*f,0.001),3.0);
      }
    }

    if(iter>=N_WARMUP){
      const R0=Math.exp(cur[0]),gamma=Math.exp(cur[1]),
            rho=Math.exp(cur[2]),phi=Math.exp(cur[3]),beta=R0*gamma;
      const It=solveSIR(beta,gamma,N,cases.length,I0);
      const mu=It.map(i=>Math.max(rho*i,1e-10));
      draws.push({R0,gamma,rho,phi,beta,mu});
    }
  }

  // Final per-parameter acceptance rates
  const rates=acc.map((a,i)=>(a/att[i]*100).toFixed(1)+'%');
  process.stdout.write(`accept[R0=${rates[0]},γ=${rates[1]},ρ=${rates[2]},φ=${rates[3]}] `);
  return draws;
}

// ═══════════════════════════════════════════════════════════════
// STATISTICS
// ═══════════════════════════════════════════════════════════════
const avg  = a=>a.reduce((s,v)=>s+v,0)/a.length;
const vari = a=>{const m=avg(a);return a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1);};

function pctile(arr,p){
  const s=[...arr].sort((a,b)=>a-b);
  const i=(s.length-1)*p;
  return s[Math.floor(i)]+(s[Math.ceil(i)]-s[Math.floor(i)])*(i-Math.floor(i));
}
function rhat(chains){
  const n=chains[0].length,m=chains.length;
  const cm=chains.map(avg),gm=avg(cm);
  const B=n/(m-1)*cm.reduce((s,c)=>s+(c-gm)**2,0);
  const W=avg(chains.map(vari));
  if(W===0)return 1;
  return Math.sqrt(((n-1)/n*W+B/n)/W);
}
function ess(chain){
  const n=chain.length,m=avg(chain),v=vari(chain);
  if(v===0)return n;
  let sum=0;
  for(let lag=1;lag<Math.min(n,300);lag++){
    let ac=0;
    for(let i=0;i<n-lag;i++) ac+=(chain[i]-m)*(chain[i+lag]-m);
    ac/=(n-lag)*v;
    if(ac<0.05)break; sum+=ac;
  }
  return Math.min(n,n/(1+2*sum));
}
function summarise(allS,key,chainArr){
  const all=allS.map(s=>s[key]);
  const chains=chainArr.map(ch=>ch.map(s=>s[key]));
  return {
    mean:avg(all), sd:Math.sqrt(vari(all)),
    q2_5:pctile(all,.025),q50:pctile(all,.5),q97_5:pctile(all,.975),
    rhat:rhat(chains), ess:ess(all)
  };
}

// ═══════════════════════════════════════════════════════════════
// PLOTLY HTML
// ═══════════════════════════════════════════════════════════════
function html(title,traces,layout={}){
  return `<!DOCTYPE html><html><head><meta charset="utf-8"><title>${title}</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>body{font-family:Georgia,serif;margin:24px;background:#fff;color:#222}
h2{margin-bottom:4px}p{color:#666;font-size:13px;margin:0 0 12px}</style>
</head><body><h2>${title}</h2>
<div id="p"></div>
<script>Plotly.newPlot('p',${JSON.stringify(traces)},
${JSON.stringify({width:900,height:480,paper_bgcolor:'#fff',
  plot_bgcolor:'#f9f9f9',font:{family:'Georgia,serif',size:12},
  legend:{bgcolor:'rgba(255,255,255,0.8)',bordercolor:'#ccc',borderwidth:1},...layout})});
</script></body></html>`;
}
function savePlot(dir,name,title,traces,layout={}){
  fs.writeFileSync(path.join(dir,name),html(title,traces,layout));
  console.log(`    ✓ ${name}`);
}

// NegBin RNG for posterior predictive
function gammaSample(shape,rate,rng){
  if(shape<1)return gammaSample(1+shape,rate,rng)*Math.pow(rng.rand(),1/shape);
  const d=shape-1/3,c=1/Math.sqrt(9*d);
  for(;;){
    let x,v; do{x=rng.randn();v=1+c*x;}while(v<=0);
    v=v**3; const u=rng.rand();
    if(u<1-0.0331*x**4)return d*v/rate;
    if(Math.log(u)<0.5*x**2+d*(1-v+Math.log(v)))return d*v/rate;
  }
}
function nbRng(mu,phi,rng){
  const lambda=gammaSample(phi,phi/Math.max(mu,1e-10),rng);
  if(lambda<=0||!isFinite(lambda))return 0;
  let L=Math.exp(-Math.min(lambda,700)),k=0,p=1;
  do{k++;p*=rng.rand();}while(p>L&&k<10000);
  return k-1;
}

// ═══════════════════════════════════════════════════════════════
// MAIN LOOP
// ═══════════════════════════════════════════════════════════════
const allResults={};

for(const site of SITES){
  console.log(`\n${'═'.repeat(56)}`);
  console.log(`  SITE: ${site.label}  —  ${site.desc}`);
  console.log(`${'═'.repeat(56)}`);

  const plotDir=path.join(PROJ,'plots',site.label);
  fs.mkdirSync(plotDir,{recursive:true});

  const rows  = readCSV(path.join(PROJ,'data',site.csv));
  const cases = rows.map(r=>Math.round(Number(r.cases)));
  const dates = rows.map(r=>r.date);
  const T=cases.length, N=site.pop, I0=site.I0;

  if(cases.every(c=>c===0)){console.log('  SKIP: all-zero');continue;}

  const peakCases=Math.max(...cases);
  const peakWeek =cases.indexOf(peakCases);
  console.log(`  Weeks:${T}  Total:${cases.reduce((a,b)=>a+b)}  Peak:${peakCases}(wk${peakWeek})  N:${N}  I0:${I0}`);

  // ── Plot 01: Raw data ──────────────────────────────────────
  savePlot(plotDir,'01_raw_data.html',
    `[${site.label}] Weekly Case Counts`,
    [{x:dates,y:cases,type:'scatter',mode:'lines+markers',
      line:{color:'#222',width:1.5},marker:{size:5,color:'#000'},name:'Cases'}],
    {xaxis:{title:'Date'},yaxis:{title:'Reported cases'}});

  // ── Plot 03: Prior predictive ──────────────────────────────
  const rngPrior=makeRng(42);
  const ppc0=[];
  for(let d=0;d<200;d++){
    const R0s=Math.exp(Math.log(1.7)+0.25*rngPrior.randn());
    const gs =Math.exp(Math.log(0.44)+0.30*rngPrior.randn());
    const rs =Math.exp(Math.log(0.05)+0.50*rngPrior.randn());
    if(rs>=1)continue;
    try{
      const It=solveSIR(R0s*gs,gs,N,T,I0);
      ppc0.push({x:dates,y:It.map(i=>rs*i),type:'scatter',mode:'lines',
        line:{color:'rgba(150,150,150,0.2)',width:0.8},showlegend:false});
    }catch{}
  }
  ppc0.push({x:dates,y:cases,type:'scatter',mode:'markers',
    marker:{color:'#000',size:6},name:'Observed'});
  savePlot(plotDir,'03_prior_predictive.html',
    `[${site.label}] Prior Predictive Check`, ppc0,
    {xaxis:{title:'Date'},yaxis:{title:'Expected cases',range:[0,peakCases*6]}});

  // ── MCMC ──────────────────────────────────────────────────
  console.log(`  Running MCMC — ${N_CHAINS} chains × ${N_ITER} iter (warmup ${N_WARMUP}):`);
  const chainSamples=[];
  for(let c=0;c<N_CHAINS;c++){
    process.stdout.write(`    Chain ${c+1}: `);
    const s=runChain(cases,N,I0,(c+1)*7919+site.label.charCodeAt(0)*31);
    chainSamples.push(s);
    console.log(`n=${s.length}`);
  }
  const allS=chainSamples.flat();
  console.log(`  Total post-warmup draws: ${allS.length}`);

  // ── Diagnostics ────────────────────────────────────────────
  console.log('\n  Convergence Diagnostics:');
  const PARAMS=['R0','beta','gamma','rho','phi'];
  const diag={};
  for(const p of PARAMS){
    diag[p]=summarise(allS,p,chainSamples);
    const d=diag[p];
    console.log(`    ${p.padEnd(6)} q50=${d.q50.toFixed(4).padStart(8)} ` +
      `95%CI=[${d.q2_5.toFixed(3)},${d.q97_5.toFixed(3)}] ` +
      `Rhat=${d.rhat.toFixed(3)} ESS=${Math.round(d.ess)}`);
  }
  const rhatOK=PARAMS.every(p=>diag[p].rhat<1.05);
  const essOK =PARAMS.every(p=>diag[p].ess >400);
  console.log(`  Convergence: R-hat ${rhatOK?'PASS ✓':'WARN ✗'}  ESS ${essOK?'PASS ✓':'WARN ✗'}`);

  // ── Plot 04: Trace plots ───────────────────────────────────
  const colours=['#222','#555','#888','#aaa'];
  const traceT=['R0','gamma','rho','phi'].flatMap((p,pi)=>
    chainSamples.map((ch,ci)=>({
      y:ch.map(s=>s[p]),x:[...Array(ch.length).keys()].map(i=>i+N_WARMUP+1),
      type:'scatter',mode:'lines',line:{color:colours[ci],width:0.7},
      name:`Chain ${ci+1}`,xaxis:`x${pi||''}`,yaxis:`y${pi||''}`,showlegend:pi===0
    }))
  );
  savePlot(plotDir,'04_trace_plots.html',
    `[${site.label}] MCMC Trace Plots`, traceT,
    {height:700,grid:{rows:4,columns:1,pattern:'independent'}});

  // ── Plot 05: Posterior histograms ──────────────────────────
  savePlot(plotDir,'05_posterior_histograms.html',
    `[${site.label}] Posterior Distributions (4000 draws × 4 chains)`,
    ['R0','gamma','rho','phi','beta'].map((p,i)=>({
      x:allS.map(s=>s[p]),type:'histogram',nbinsx:60,name:p,
      marker:{color:'#777',line:{color:'#333',width:0.4}},
      xaxis:`x${i||''}`,yaxis:`y${i||''}`,showlegend:false
    })),
    {height:600,grid:{rows:2,columns:3,pattern:'independent'}});

  // ── Posterior predictive ───────────────────────────────────
  const muMat=allS.map(s=>s.mu);
  const muMed =Array.from({length:T},(_,t)=>pctile(muMat.map(r=>r[t]),.50));
  const muLo90=Array.from({length:T},(_,t)=>pctile(muMat.map(r=>r[t]),.05));
  const muHi90=Array.from({length:T},(_,t)=>pctile(muMat.map(r=>r[t]),.95));
  const muLo50=Array.from({length:T},(_,t)=>pctile(muMat.map(r=>r[t]),.25));
  const muHi50=Array.from({length:T},(_,t)=>pctile(muMat.map(r=>r[t]),.75));

  const casesAvg=avg(cases);
  const ssTot=cases.reduce((s,c)=>s+(c-casesAvg)**2,0);
  const ssRes=cases.reduce((s,c,t)=>s+(c-muMed[t])**2,0);
  const R2=ssTot>0?1-ssRes/ssTot:0;
  const RMSE=Math.sqrt(ssRes/T);
  console.log(`  R² = ${R2.toFixed(3)}  RMSE = ${RMSE.toFixed(2)}`);

  // ── Plot 07: PPC ──────────────────────────────────────────
  savePlot(plotDir,'07_posterior_predictive.html',
    `[${site.label}] Posterior Predictive Check  (R²=${R2.toFixed(2)})`,
    [
      {x:[...dates,...[...dates].reverse()],
       y:[...muHi90,...[...muLo90].reverse()],
       fill:'toself',fillcolor:'rgba(190,190,190,0.4)',
       line:{color:'transparent'},name:'90% CI'},
      {x:[...dates,...[...dates].reverse()],
       y:[...muHi50,...[...muLo50].reverse()],
       fill:'toself',fillcolor:'rgba(110,110,110,0.5)',
       line:{color:'transparent'},name:'50% CI'},
      {x:dates,y:muMed,type:'scatter',mode:'lines',
       line:{color:'#111',width:2.5},name:'Median fit'},
      {x:dates,y:cases,type:'scatter',mode:'markers',
       marker:{color:'#000',size:7,symbol:'circle-open',line:{width:2}},name:'Observed'}
    ],
    {xaxis:{title:'Date'},yaxis:{title:'Cases'}});

  // ── Plot 08: Simulated observations ───────────────────────
  const rngR=makeRng(999);
  const repT=[];
  const nRep=Math.min(120,allS.length);
  for(let d=0;d<nRep;d++){
    const ri=Math.floor(rngR.rand()*allS.length);
    repT.push({x:dates,y:allS[ri].mu.map(mu=>nbRng(mu,allS[ri].phi,rngR)),
      type:'scatter',mode:'lines',line:{color:'rgba(150,150,150,0.18)',width:0.7},showlegend:false});
  }
  repT.push({x:dates,y:cases,type:'scatter',mode:'markers',
    marker:{color:'#000',size:6,symbol:'circle-open',line:{width:2}},name:'Observed'});
  savePlot(plotDir,'08_posterior_predictive_sim.html',
    `[${site.label}] Posterior Predictive — Simulated Observations`, repT,
    {xaxis:{title:'Date'},yaxis:{title:'Cases'}});

  // ── Plot 09: R0 posterior ─────────────────────────────────
  const R0arr=allS.map(s=>s.R0), R0med=pctile(R0arr,.5);
  savePlot(plotDir,'09_R0_distribution.html',
    `[${site.label}] Posterior R₀  (median=${R0med.toFixed(2)}, 95%CI=[${pctile(R0arr,.025).toFixed(2)},${pctile(R0arr,.975).toFixed(2)}])`,
    [
      {x:R0arr,type:'histogram',nbinsx:60,
       marker:{color:'#888',line:{color:'#333',width:0.5}},name:'Posterior'},
      {x:[1,1],y:[0,allS.length*.12],type:'scatter',mode:'lines',
       line:{color:'#333',width:2,dash:'dot'},name:'R₀ = 1 threshold'},
      {x:[R0med,R0med],y:[0,allS.length*.12],type:'scatter',mode:'lines',
       line:{color:'#111',width:2},name:'Median'}
    ],
    {xaxis:{title:'R₀'},yaxis:{title:'Frequency'}});

  // ── Plot 10: Residuals ────────────────────────────────────
  const resid=cases.map((c,t)=>c-muMed[t]);
  savePlot(plotDir,'10_residuals.html',
    `[${site.label}] Residuals (Observed − Fitted)`,
    [
      {x:dates,y:Array(T).fill(0),type:'scatter',mode:'lines',
       line:{color:'#bbb',dash:'dash'},showlegend:false},
      {x:dates,y:resid,type:'scatter',mode:'lines+markers',
       line:{color:'#444',width:1.2},marker:{size:5,color:'#222'},name:'Residual'}
    ],
    {xaxis:{title:'Date'},yaxis:{title:'Observed − Fitted'}});

  // ── Plot 06: β-γ scatter ──────────────────────────────────
  savePlot(plotDir,'06_beta_gamma_scatter.html',
    `[${site.label}] Posterior β vs γ  (Pearson r=${
      (()=>{const bv=allS.map(s=>s.beta),gv=allS.map(s=>s.gamma),
            mb=avg(bv),mg=avg(gv);
            return (allS.reduce((s,_,i)=>s+(bv[i]-mb)*(gv[i]-mg),0)/allS.length
              /Math.sqrt(vari(bv)*vari(gv))).toFixed(3);})()})`,
    [{x:allS.map(s=>s.gamma),y:allS.map(s=>s.beta),type:'scatter',mode:'markers',
      marker:{size:2,color:'rgba(60,60,60,0.25)'},name:'Posterior draws'}],
    {xaxis:{title:'γ (recovery rate)'},yaxis:{title:'β (transmission rate)'}});

  // ── Plot 02: SIR compartments with posterior median ───────
  const bM=diag.beta.q50,gM=diag.gamma.q50,rM=diag.rho.q50;
  const ItM=solveSIR(bM,gM,N,T,I0);
  savePlot(plotDir,'02_sir_compartments.html',
    `[${site.label}] SIR Compartments (posterior medians: R₀=${diag.R0.q50.toFixed(2)})`,
    [
      {x:dates,y:ItM.map(i=>(N-I0)-ItM.reduce((s,v,j)=>j<=dates.indexOf(dates[0])?s:s,0)),
       type:'scatter',mode:'lines',line:{color:'#222',width:1.5,dash:'solid'},name:'Susceptible (scaled)'},
      {x:dates,y:ItM,type:'scatter',mode:'lines',
       line:{color:'#555',width:1.5,dash:'dash'},name:'Infectious I(t)'},
      {x:dates,y:ItM.map(i=>rM*i),type:'scatter',mode:'lines',
       line:{color:'#888',width:1.5,dash:'dot'},name:'Expected cases ρ·I(t)'},
      {x:dates,y:cases,type:'scatter',mode:'markers',
       marker:{color:'#000',size:5,symbol:'circle-open',line:{width:2}},name:'Observed'}
    ],
    {xaxis:{title:'Date'},yaxis:{title:'Count / Cases'}});

  // ── Save JSON ─────────────────────────────────────────────
  const result={
    site:site.label, description:site.desc, N, I0, n_weeks:T,
    total_cases:cases.reduce((a,b)=>a+b), peak_cases:peakCases, peak_week:peakWeek,
    n_draws:allS.length, convergence:{rhat_pass:rhatOK,ess_pass:essOK},
    parameters:Object.fromEntries(PARAMS.map(p=>[p,{
      mean:+diag[p].mean.toFixed(4), sd:+diag[p].sd.toFixed(4),
      q2_5:+diag[p].q2_5.toFixed(4), q50:+diag[p].q50.toFixed(4),
      q97_5:+diag[p].q97_5.toFixed(4),
      rhat:+diag[p].rhat.toFixed(4), ess:Math.round(diag[p].ess)
    }])),
    fit:{R2:+R2.toFixed(4),RMSE:+RMSE.toFixed(3)}
  };
  fs.writeFileSync(path.join(OUT,`results_${site.label}.json`),JSON.stringify(result,null,2));
  console.log(`  ✓ outputs/results_${site.label}.json`);
  console.log(`  ✓ plots/${site.label}/  (10 HTML files)`);
  allResults[site.label]=result;
}

// ═══════════════════════════════════════════════════════════════
// CROSS-SITE SUMMARY
// ═══════════════════════════════════════════════════════════════
console.log(`\n\n${'═'.repeat(64)}`);
console.log('  CROSS-SITE COMPARISON');
console.log(`${'═'.repeat(64)}`);
console.log('  Site    R₀ median  95% CI          β       γ       ρ      R²');
console.log('  '+'─'.repeat(60));
for(const[nm,r] of Object.entries(allResults)){
  const p=r.parameters;
  console.log(
    `  ${nm.padEnd(7)} ${p.R0.q50.toFixed(3).padStart(8)}  `+
    `[${p.R0.q2_5.toFixed(2)},${p.R0.q97_5.toFixed(2)}]`.padEnd(16)+
    ` ${p.beta.q50.toFixed(3).padStart(6)}  ${p.gamma.q50.toFixed(3).padStart(6)} `+
    ` ${p.rho.q50.toFixed(4).padStart(7)}  ${r.fit.R2.toFixed(3).padStart(5)}`
  );
}

fs.writeFileSync(path.join(OUT,'results_all_sites.json'),JSON.stringify(allResults,null,2));
console.log('\n  ✓ outputs/results_all_sites.json');
console.log('\n=== All sites complete ===\n');
