(* ::Package:: *)

(*
  1-loop massless cutbubble at eps = 1/2 (D = 3) — Phase 3 oracle
  for the ε-extremes axis (boundary case at large eps).  Same
  cutbubble family as cutbubble_1L_eps001 but evaluated at
  eps = 1/2 (D = 4 - 2*eps = 3 dimensions).  Note: an initial
  attempt at eps = 1 (D = 2) hit an upstream MMA AMFlow limitation
  where the DESolver returned a partially-symbolic expression
  `(1/2π) Im[DESolver`Private`variables[1, 1]]`, so eps = 1/2 was
  chosen as the next-most-extreme rational eps that yields a clean
  numeric result.  Tests the Laurent-expansion machinery at the
  far end of dimensional regularization where 1/eps^n coefficients
  are O(1) rather than asymptotically large.  Target is
  j[cutbubble, 1, 1] at s = 4, eps = 1/2.
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions[
  "IBPReducer" -> "Kira",
  "BlackBoxRank" -> 0,
  "BlackBoxDot" -> 1,
  "DeleteBlackBoxDirectory" -> False
];
Kira`$KiraExecutable = "/usr/local/bin/kira";
Kira`$FermatExecutable = "/usr/share/Ferl7/fer64";

SetAMFOptions[
  "AMFMode" -> {"Prescription", "Mass", "Propagator"},
  "EndingScheme" -> {"Tradition", "Cutkosky", "SingleMass"},
  "WorkingPre" -> 120,
  "ChopPre" -> 20,
  "RationalizePre" -> 100,
  "XOrder" -> 200,
  "ExtraXOrder" -> 240,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = cutbubble;
AMFlowInfo["Loop"] = {l};
AMFlowInfo["Leg"] = {p};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p^2 -> s};
AMFlowInfo["Propagator"] = {
  l^2,
  (l + p)^2
};
AMFlowInfo["Cut"] = {1, 1};
AMFlowInfo["Numeric"] = {s -> 4};
AMFlowInfo["NThread"] = 4;

targets = {
  j[cutbubble, 1, 1]
};
epslist = {1/2};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "cutbubble_1L_eps2_mma_cache"}]];

Print["==BENCH== mma 1-loop cutbubble BlackBoxAMFlow eps=1/2 (D=3)"];
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
