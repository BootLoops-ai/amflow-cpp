(* ::Package:: *)

(*
  1-loop pentagon with TWO different internal masses placed on
  props 0 and 2 (mAsq on l^2, mBsq on (l+p1+p2)^2).  Completes the
  1L 5-leg mass axis: massless (pentagon_1L_eps001), 1-mass
  (pentagon_1L_W_mass), 2-mass (this), 3-mass (pentagon_1L_3mass).
  Target: corner j[pentagon2m, 1,1,1,1,1] at eps = 1/1000.
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

AMFlowInfo["Family"] = pentagon2m;
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
  (l + p1 + p2 + p3 + p4)^2
};
AMFlowInfo["Numeric"] = {
  s12 -> -2, s23 -> -3, s34 -> -5,
  s13 -> -2, s14 -> 15, s24 -> -3,
  mAsq -> 1, mBsq -> 4
};
AMFlowInfo["NThread"] = 4;

targets = { j[pentagon2m, 1, 1, 1, 1, 1] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "pentagon_1L_2mass_eps001_mma_cache"}]];

Print["==BENCH== mma 1-loop pentagon with 2 distinct masses"];
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
