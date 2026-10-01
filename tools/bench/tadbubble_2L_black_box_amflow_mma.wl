(* ::Package:: *)

(*
  Factorizable 2-loop topology: 1-loop bubble carrying the external
  momentum (l1) and a 1-loop tadpole (l2) hanging on a separate loop
  variable.  The graph factors as (bubble)*(tadpole), so the IBP
  reduction and AMFlow path see a degenerate sub-system.  Distinctive:
  separable / factorizable 2L family.
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

AMFlowInfo["Family"] = tadbub;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p^2 -> s};
AMFlowInfo["Propagator"] = {
  l1^2 - msq,
  (l1 + p)^2 - msq,
  l2^2 - msq,
  (l1 + l2)^2,
  (l2 + p)^2
};
AMFlowInfo["Numeric"] = {s -> 4, msq -> 1};
AMFlowInfo["NThread"] = 4;

targets = {j[tadbub, 1, 1, 1, 0, 0]};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "tadbubble_2L_mma_cache"}]];

Print["==BENCH== mma 2-loop factorizable bubble+tadpole BlackBoxAMFlow"];
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
