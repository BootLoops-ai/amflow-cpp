(* ::Package:: *)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions[
  "IBPReducer" -> "Kira",
  "BlackBoxRank" -> 3,
  "BlackBoxDot" -> 0,
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

AMFlowInfo["Family"] = tt;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1 - p2 - p3};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0,
  p2^2 -> 0,
  p3^2 -> msq,
  p4^2 -> msq,
  (p1 + p2)^2 -> s,
  (p1 + p3)^2 -> t
};
AMFlowInfo["Propagator"] = {
  l1^2,
  (l1 + p1)^2,
  (l1 + p1 + p2)^2,
  l2^2,
  -msq + (l2 + p3)^2,
  (l2 + p3 + p4)^2,
  (l1 + l2)^2,
  (l1 - p3)^2,
  (l2 + p1)^2
};
AMFlowInfo["Numeric"] = {s -> 30, t -> -10/3, msq -> 1};
AMFlowInfo["NThread"] = 4;

targets = {
  j[tt, 1, 1, 1, 1, 1, 1, 1, -3, 0],
  j[tt, 1, 1, 1, 1, 1, 1, 1, -2, -1],
  j[tt, 1, 1, 1, 1, 1, 1, 1, -1, -2],
  j[tt, 1, 1, 1, 1, 1, 1, 1, 0, -3]
};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = "/tmp/desolver_bench_tt_2loop_box_mma"];

Print["==BENCH== mma tt 2-loop box BlackBoxAMFlow sampled"];
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
