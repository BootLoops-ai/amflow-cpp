(* ::Package:: *)

(*
  2L planar vertex 4mass at eps=1/100.
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

AMFlowInfo["Family"] = vtx24me01;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0,
  p2^2 -> 0,
  (p1 + p2)^2 -> s
};
AMFlowInfo["Propagator"] = {
  l1^2 - mAsq,
  (l1 + p1)^2 - mBsq,
  (l1 + p1 + p2)^2 - mCsq,
  l2^2 - mDsq,
  (l2 + p1 + p2)^2,
  (l1 - l2)^2,
  (l2 + p1)^2
};
AMFlowInfo["Numeric"] = {
  s -> -1,
  mAsq -> 1,
  mBsq -> 4,
  mCsq -> 9,
  mDsq -> 16
};
AMFlowInfo["NThread"] = 4;

targets = { j[vtx24me01, 1, 1, 1, 1, 1, 1, 0] };
epslist = {1/100};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "vtx2_2L_4mass_eps01_mma_cache"}]];

Print["==BENCH== mma 2L planar vertex 4mass at eps=1/100."];
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
