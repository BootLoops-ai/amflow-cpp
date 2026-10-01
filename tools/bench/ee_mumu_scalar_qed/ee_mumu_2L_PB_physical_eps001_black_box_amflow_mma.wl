(* ::Package:: *)

(*
  e^- e^+ -> mu^- mu^+ scalar QED — 2-loop PLANAR double-box (PB)
  Same family structure as Bhabha PB (Henn-Smirnov-Smirnov) but with
  TWO different masses: top rail = electron (mesq), bottom rail =
  muon (mmsq).  Physical kinematics: s above 4 mmsq threshold.

  Kinematics (physical):
    p1, p2 = electrons in,  p3, p4 = muons in (all-incoming)
    p1^2 = p2^2 = mesq = 1
    p3^2 = p4^2 = mmsq = 4   (so m_mu = 2 m_e)
    s = (p1+p2)^2 = 25  (above 4 mmsq = 16 threshold)
    t = (p1+p3)^2 = -5  (in physical range [-14.36, -0.64])
    u = 10 - 25 - (-5) = -10
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions[
  "IBPReducer" -> "Kira",
  "BlackBoxRank" -> 0,
  "BlackBoxDot" -> 3,
  "DeleteBlackBoxDirectory" -> False
];
Kira`$KiraExecutable = "/usr/local/bin/kira";
Kira`$FermatExecutable = "/usr/share/Ferl7/fer64";

SetAMFOptions[
  "AMFMode" -> {"Prescription", "Mass", "Propagator"},
  "EndingScheme" -> {"Tradition", "SingleMass"},
  "WorkingPre" -> 120, "ChopPre" -> 20, "RationalizePre" -> 100,
  "XOrder" -> 240, "ExtraXOrder" -> 280,
  "LearnXOrder" -> -1, "TestXOrder" -> 5
];

AMFlowInfo["Family"] = mupairpb;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1 - p2 - p3};
AMFlowInfo["Replacement"] = {
  p1^2 -> mesq,
  p2^2 -> mesq,
  p3^2 -> mmsq,
  p4^2 -> mmsq,
  (p1 + p2)^2 -> s,
  (p1 + p3)^2 -> t
};
AMFlowInfo["Propagator"] = {
  l1^2 - mesq,
  (l1 + p1 + p2)^2 - mesq,
  l2^2 - mmsq,
  (l2 + p1 + p2)^2 - mmsq,
  (l1 + p1)^2,
  (l1 - l2)^2,
  (l2 + p3)^2,
  (l2 + p1)^2,
  (l1 + p3)^2
};
AMFlowInfo["Numeric"] = {
  mesq -> 1,
  mmsq -> 4,
  s -> 25,
  t -> -5
};
AMFlowInfo["NThread"] = 4;

targets = { j[mupairpb, 1, 1, 1, 1, 1, 1, 1, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "ee_mumu_scalar_qed", "mma_refs",
    "ee_mumu_2L_PB_physical_eps001_mma_cache"}]];

Print["==BENCH== mma ee->mumu 2L PB physical (mesq=1, mmsq=4, s=25, t=-5)."];
Print["==CACHE== ", cache];
Print["==EPS== ", epslist // InputForm];

elapsed = AbsoluteTiming[
  sol = BlackBoxAMFlow[targets, epslist, cache];
][[1]];

Print["==TIME== ", elapsed];
Print["==RESULT_INPUTFORM== ", sol // InputForm];

Do[
  key = Keys[sol][[i]];
  val = Values[sol][[i, 1]];
  Print["==VALUE== ", key // InputForm, " :: ",
    N[Re[val], 50] // InputForm, " :: ", N[Im[val], 50] // InputForm],
  {i, Length[Keys[sol]]}
];

sysSol = Get[FileNameJoin[{cache, "1", "solution"}]];
If[Head[sysSol] === List,
  Do[
    key = sysSol[[i, 1]];
    val = sysSol[[i, 2, 1]];
    Print["==SAMPLED== ", key // InputForm, " :: ",
      N[Re[val], 50] // InputForm, " :: ", N[Im[val], 50] // InputForm],
    {i, Length[sysSol]}
  ]];

Quit[];
