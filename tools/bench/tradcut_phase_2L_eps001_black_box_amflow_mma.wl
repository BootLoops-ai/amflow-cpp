(* ::Package:: *)

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
  "EndingScheme" -> {"Tradition", "Cutkosky", "SingleMass"},
  "WorkingPre" -> 100,
  "ChopPre" -> 20,
  "RationalizePre" -> 100,
  "XOrder" -> 200,
  "ExtraXOrder" -> 240,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

(* Tradition-with-cut benchmark: 2-loop family from upstream
   automatic_phasespace example.  Cut on 3 of 7 propagators.  Cutkosky
   does not apply (the cut topology is not phase_volume); upstream
   AMFlow runs the Tradition scheme, which is what exercises the C++
   port's `build_boundary` cut-projection (Phase 1B fix). *)

AMFlowInfo["Family"] = phase;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p1^2 -> 0, p2^2 -> 0, (p1+p2)^2 -> s};
AMFlowInfo["Propagator"] = {
  l1^2 - msq,
  (l1+p1)^2,
  l2^2,
  (l1+l2+p1)^2,
  (l1+l2+p1+p2)^2,
  (l1+l2+p2)^2,
  (l1+p2)^2
};
AMFlowInfo["Prescription"] = {0, 0};
AMFlowInfo["Cut"] = {1, 0, 1, 0, 1, 0, 0};
AMFlowInfo["Numeric"] = {s -> 100, msq -> 1};
AMFlowInfo["NThread"] = 4;

targets = {
  j[phase, 1, 0, 1, 0, 1, 0, 0]
};
epslist = {1/100};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "tradcut_phase_2L_eps001_mma_cache"}]];

Print["==BENCH== mma 2-loop Tradition-with-cut phase-space top BlackBoxAMFlow sampled"];
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
