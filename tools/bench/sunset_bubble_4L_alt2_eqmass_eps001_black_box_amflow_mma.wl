(* ::Package:: *)

(*
  4L 2-leg: bubble decoration on l3 instead of l1.
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

AMFlowInfo["Family"] = sb4alt2eq;
AMFlowInfo["Loop"] = {l1, l2, l3, l4};
AMFlowInfo["Leg"] = {p1, p2};
AMFlowInfo["Conservation"] = {p2 -> -p1};
AMFlowInfo["Replacement"] = {
  p1^2 -> psq
};
AMFlowInfo["Propagator"] = {
  l1^2 - msq,
  l2^2 - msq,
  l3^2 - msq,
  (l3 - l4)^2 - msq,
  (l1 + l2 + l3 - p1)^2 - msq,
  l4^2 - msq,
  (l1 + p1)^2,
  (l2 + p1)^2,
  (l3 + p1)^2,
  (l4 + p1)^2,
  (l1 - l2)^2,
  (l1 - l3)^2,
  (l1 - l4)^2,
  (l2 - l4)^2
};
AMFlowInfo["Numeric"] = {
  psq -> 1,
  msq -> 1
};
AMFlowInfo["NThread"] = 4;

targets = { j[sb4alt2eq, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "sunset_bubble_4L_alt2_eqmass_eps001_mma_cache"}]];

Print["==BENCH== mma 4L 2-leg: bubble decoration on l3 instead of l1."];
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
