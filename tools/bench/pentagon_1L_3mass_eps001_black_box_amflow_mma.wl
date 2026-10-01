(* ::Package:: *)

(*
  1-loop pentagon with THREE different internal masses placed on
  non-cyclically-adjacent propagators (props 0, 2, 4 of the 5-prop
  basis): mAsq on l^2, mBsq on (l+p1+p2)^2, mCsq on
  (l+p1+p2+p3+p4)^2.  Pushes the 3-mass-axis at 1L 5-leg —
  previously the suite has 1L 5-leg only at single-mass (pentagon_1L_W_mass)
  or all-massless (pentagon_1L_eps001).

  Massive propagators carry distinct hierarchical masses
  (mAsq=1, mBsq=4, mCsq=9) so the SingleMass / Mass scheme dispatch
  has to discriminate three mass values along the η-injection axis.

  Same kinematic numerics as pentagon_1L_eps001 for direct
  comparability against the massless oracle.
  Target: corner j[pentagon3m, 1,1,1,1,1] at eps = 1/1000.
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
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

AMFlowInfo["Family"] = pentagon3m;
AMFlowInfo["Loop"] = {l};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
  (p1 + p2)^2 -> s12, (p1 + p3)^2 -> s13, (p1 + p4)^2 -> s14,
  (p2 + p3)^2 -> s23, (p2 + p4)^2 -> s24, (p3 + p4)^2 -> s34
};
AMFlowInfo["Propagator"] = {
  l^2 - mAsq,
  (l + p1)^2,
  (l + p1 + p2)^2 - mBsq,
  (l + p1 + p2 + p3)^2,
  (l + p1 + p2 + p3 + p4)^2 - mCsq
};
AMFlowInfo["Numeric"] = {
  s12 -> -2, s23 -> -3, s34 -> -5,
  s13 -> -2, s14 -> 15, s24 -> -3,
  mAsq -> 1, mBsq -> 4, mCsq -> 9
};
AMFlowInfo["NThread"] = 4;

targets = { j[pentagon3m, 1, 1, 1, 1, 1] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "pentagon_1L_3mass_eps001_mma_cache"}]];

Print["==BENCH== mma 1-loop pentagon with 3 distinct masses"];
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

Quit[];
