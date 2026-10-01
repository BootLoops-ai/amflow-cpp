(* ::Package:: *)

(*
  1-loop pentagon — Phase 3 oracle (external-leg axis = 5).
  All-massless internal lines, 5 external on-shell legs (p1..p5 with
  p5 = -(p1+p2+p3+p4) via conservation), 5 cyclic propagators.
  Kinematic basis: s12=-2, s23=-3, s34=-5, s45=-7, s51=-11 (with
  derived s13=-2, s14=+15, s24=-3, s25=+8, s35=+10 from 5-pt
  Mandelstam closure).  Target: corner j[pentagon, 1, 1, 1, 1, 1]
  at eps = 1/1000.
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
  "WorkingPre" -> 120,
  "ChopPre" -> 20,
  "RationalizePre" -> 100,
  "XOrder" -> 240,
  "ExtraXOrder" -> 280,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = pentagon;
AMFlowInfo["Loop"] = {l};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
  (p1 + p2)^2 -> s12,
  (p1 + p3)^2 -> s13,
  (p1 + p4)^2 -> s14,
  (p2 + p3)^2 -> s23,
  (p2 + p4)^2 -> s24,
  (p3 + p4)^2 -> s34
};
AMFlowInfo["Propagator"] = {
  l^2,
  (l + p1)^2,
  (l + p1 + p2)^2,
  (l + p1 + p2 + p3)^2,
  (l + p1 + p2 + p3 + p4)^2
};
AMFlowInfo["Numeric"] = {
  s12 -> -2, s23 -> -3, s34 -> -5,
  s13 -> -2, s14 -> 15, s24 -> -3
};
AMFlowInfo["NThread"] = 4;

targets = { j[pentagon, 1, 1, 1, 1, 1] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "pentagon_1L_eps001_mma_cache"}]];

Print["==BENCH== mma 1-loop massless pentagon BlackBoxAMFlow sampled"];
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
