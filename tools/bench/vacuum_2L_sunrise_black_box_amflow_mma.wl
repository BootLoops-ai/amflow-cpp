(* ::Package:: *)

(*
  2-loop pure vacuum integral (NO external legs).  Equal-mass sunrise
  topology.  Exercises the AMFlow Vacuum lookup table path (auto-applied
  Vacuum[2, 3] table entry) rather than the standard differential-equation
  pipeline.  Distinctive: 0 external momenta.
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

AMFlowInfo["Family"] = vac2;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {};
AMFlowInfo["Propagator"] = {
  l1^2 - msq,
  l2^2 - msq,
  (l1 + l2)^2 - msq
};
AMFlowInfo["Numeric"] = {msq -> 1};
AMFlowInfo["NThread"] = 4;

targets = {j[vac2, 1, 1, 1]};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "vacuum_2L_sunrise_mma_cache"}]];

Print["==BENCH== mma 2-loop pure vacuum sunrise BlackBoxAMFlow"];
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
