(* ::Package:: *)

(*
  Cut 4L banana with 2 masses.
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

AMFlowInfo["Family"] = cut42m;
AMFlowInfo["Loop"] = {l1, l2, l3, l4};
AMFlowInfo["Leg"] = {p};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {
  p^2 -> s
};
AMFlowInfo["Propagator"] = {
  l1^2 - mAsq,
  l2^2 - mBsq,
  l3^2,
  l4^2,
  (l1 + l2 + l3 + l4 - p)^2,
  (l1 + p)^2,
  (l2 + p)^2,
  (l3 + p)^2,
  (l4 + p)^2,
  (l1 - l2)^2,
  (l1 - l3)^2,
  (l1 - l4)^2,
  (l2 - l3)^2,
  (l2 - l4)^2
};
AMFlowInfo["Numeric"] = {
  s -> 10,
  mAsq -> 1,
  mBsq -> 4
};
AMFlowInfo["NThread"] = 4;

targets = { j[cut42m, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "cutbanana_4L_2mass_eps001_mma_cache"}]];

Print["==BENCH== mma Cut 4L banana with 2 masses."];
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
