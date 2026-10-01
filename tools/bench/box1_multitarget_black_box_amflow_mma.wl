(* ::Package:: *)

(*
  1-loop massless box with THREE targets in a single BlackBoxAMFlow
  call: corner + two dotted variants.  Distinctive: multi-target
  output in one call (every prior bench had one primary target).
  Exercises the "for each requested target" loop in the C++
  amflow_cli output formatter.
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
  "EndingScheme" -> {"Tradition", "SingleMass"},
  "WorkingPre" -> 120,
  "ChopPre" -> 20,
  "RationalizePre" -> 100,
  "XOrder" -> 200,
  "ExtraXOrder" -> 240,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = box1mt;
AMFlowInfo["Loop"] = {l};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1 - p2 - p3};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
  (p1 + p2)^2 -> s, (p1 + p3)^2 -> t
};
AMFlowInfo["Propagator"] = {
  l^2,
  (l + p1)^2,
  (l + p1 + p2)^2,
  (l + p1 + p2 + p4)^2
};
AMFlowInfo["Numeric"] = {s -> -3, t -> -1};
AMFlowInfo["NThread"] = 4;

targets = {
  j[box1mt, 1, 1, 1, 1],
  j[box1mt, 2, 1, 1, 1],
  j[box1mt, 1, 1, 2, 1]
};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "box1_multitarget_mma_cache"}]];

Print["==BENCH== mma 1-loop box multi-target BlackBoxAMFlow"];
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
    N[Re[val], 40] // InputForm, " :: ", N[Im[val], 40] // InputForm],
  {i, Length[Keys[sol]]}
];

Quit[];
