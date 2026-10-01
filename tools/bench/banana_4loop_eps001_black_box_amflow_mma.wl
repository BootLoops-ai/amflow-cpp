(* ::Package:: *)

(*
  4-loop equal-mass banana sunrise — Phase 3 oracle (loop-number
  axis L=4).  Family: 5 massive propagators (each loop momentum
  squared minus msq, plus a "wrapping" propagator gluing the four
  loops back to p1) and 9 ISPs spanning the remaining
  loop-momenta scalar products.  Target is the equal-mass top
  sector j[banana4, 1, 1, 1, 1, 1, 0, ..., 0] at psq = -3, msq = 1,
  eps = 1/1000.
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions[
  "IBPReducer" -> "Kira",
  "BlackBoxRank" -> 0,
  "BlackBoxDot" -> 5,
  "DeleteBlackBoxDirectory" -> False
];
Kira`$KiraExecutable = "/usr/local/bin/kira";
Kira`$FermatExecutable = "/usr/share/Ferl7/fer64";

SetAMFOptions[
  "AMFMode" -> {"Prescription", "Mass", "Propagator"},
  "EndingScheme" -> {"Tradition", "SingleMass"},
  "WorkingPre" -> 120,
  "ChopPre" -> 20,
  "RationalizePre" -> 100,
  "XOrder" -> 240,
  "ExtraXOrder" -> 280,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = banana4;
AMFlowInfo["Loop"] = {l1, l2, l3, l4};
AMFlowInfo["Leg"] = {p1};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p1^2 -> psq};
AMFlowInfo["Propagator"] = {
  l1^2 - msq,
  l2^2 - msq,
  l3^2 - msq,
  l4^2 - msq,
  (l1 + l2 + l3 + l4 - p1)^2 - msq,
  (l1 + p1)^2,
  (l2 + p1)^2,
  (l3 + p1)^2,
  (l4 + p1)^2,
  (l1 - l2)^2,
  (l1 - l3)^2,
  (l1 - l4)^2,
  (l2 - l3)^2,
  (l2 - l4)^2
};
AMFlowInfo["Numeric"] = {psq -> -3, msq -> 1};
AMFlowInfo["NThread"] = 4;

targets = {
  j[banana4, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "banana_4loop_eps001_mma_cache"}]];

Print["==BENCH== mma 4-loop equal-mass banana BlackBoxAMFlow sampled"];
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
