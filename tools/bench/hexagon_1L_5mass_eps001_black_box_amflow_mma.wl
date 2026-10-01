(* ::Package:: *)

(*
  1-loop hexagon with 5 distinct masses on props 0..4 (5 of 6).
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

AMFlowInfo["Family"] = hex15m;
AMFlowInfo["Loop"] = {l};
AMFlowInfo["Leg"] = {p1, p2, p3, p4, p5};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0,
  p2^2 -> 0,
  p3^2 -> 0,
  p4^2 -> 0,
  p5^2 -> 0,
  (p1 + p2)^2 -> s12,
  (p1 + p3)^2 -> s13,
  (p1 + p4)^2 -> s14,
  (p1 + p5)^2 -> s15,
  (p2 + p3)^2 -> s23,
  (p2 + p4)^2 -> s24,
  (p2 + p5)^2 -> s25,
  (p3 + p4)^2 -> s34,
  (p3 + p5)^2 -> s35,
  (p4 + p5)^2 -> s45
};
AMFlowInfo["Propagator"] = {
  l^2 - mAsq,
  (l + p1)^2 - mBsq,
  (l + p1 + p2)^2 - mCsq,
  (l + p1 + p2 + p3)^2 - mDsq,
  (l + p1 + p2 + p3 + p4)^2 - mEsq,
  (l + p1 + p2 + p3 + p4 + p5)^2
};
AMFlowInfo["Numeric"] = {
  s12 -> -2,
  s13 -> 1,
  s14 -> 3,
  s15 -> 5,
  s23 -> -3,
  s24 -> -2,
  s25 -> 1,
  s34 -> -5,
  s35 -> 4,
  s45 -> -2,
  mAsq -> 1,
  mBsq -> 4,
  mCsq -> 9,
  mDsq -> 16,
  mEsq -> 25
};
AMFlowInfo["NThread"] = 4;

targets = { j[hex15m, 1, 1, 1, 1, 1, 1] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "hexagon_1L_5mass_eps001_mma_cache"}]];

Print["==BENCH== mma 1-loop hexagon with 5 distinct masses on props 0..4 (5 of 6)."];
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
