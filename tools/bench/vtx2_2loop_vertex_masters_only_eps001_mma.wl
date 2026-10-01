(* ::Package:: *)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions[
  "IBPReducer" -> "Kira",
  "BlackBoxRank" -> 2,
  "BlackBoxDot" -> 0,
  "DeleteBlackBoxDirectory" -> False
];
Kira`$KiraExecutable = "/usr/local/bin/kira";
Kira`$FermatExecutable = "/usr/share/Ferl7/fer64";

SetAMFOptions[
  "AMFMode" -> {"Prescription", "Mass", "Propagator"},
  "EndingScheme" -> {"Tradition", "SingleMass"},
  "WorkingPre" -> 100,
  "ChopPre" -> 20,
  "RationalizePre" -> 80,
  "XOrder" -> 100,
  "ExtraXOrder" -> 120,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = vtx2;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0,
  p2^2 -> 0,
  (p1 + p2)^2 -> s
};
AMFlowInfo["Propagator"] = {
  l1^2,
  (l1 + p1)^2,
  (l1 + p1 + p2)^2,
  l2^2,
  (l2 + p1 + p2)^2,
  (l1 - l2)^2,
  (l2 + p1)^2
};
AMFlowInfo["Numeric"] = {s -> -1};
AMFlowInfo["NThread"] = 4;

targets = {
  j[vtx2, 0, 0, 1, 1, 0, 1, 0],
  j[vtx2, 1, 0, 1, 1, 1, 0, 0],
  j[vtx2, 0, 1, 0, 1, 1, 1, 0]
};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = "/tmp/desolver_bench_vtx2_masters_only_mma"];

Print["==BENCH== mma vtx2 masters-only BlackBoxAMFlow sampled"];
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
