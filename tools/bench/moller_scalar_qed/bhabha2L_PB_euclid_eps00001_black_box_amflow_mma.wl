(* ::Package:: *)

(*
  Bhabha 2-loop planar double-box (PB) for massive Bhabha scattering.

  Topology: arXiv:1307.4083 (Henn, Smirnov, Smirnov, JHEP 11(2013)041) eq. (2.1).

  7-propagator parent topology; AMFlow padded with 2 ISPs to rank 9.
  Region: Euclidean-like (s, t both spacelike; u in u-channel physical region).
  eps point: 1/1000.
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", "..", ".."}];
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
  "WorkingPre" -> 200, "ChopPre" -> 20, "RationalizePre" -> 100,
  "XOrder" -> 400, "ExtraXOrder" -> 480,
  "LearnXOrder" -> -1, "TestXOrder" -> 5
];

AMFlowInfo["Family"] = bhabha2pb;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1 - p2 - p3};
AMFlowInfo["Replacement"] = {
  p1^2 -> msq,
  p2^2 -> msq,
  p3^2 -> msq,
  p4^2 -> msq,
  (p1 + p2)^2 -> s,
  (p1 + p3)^2 -> t
};
AMFlowInfo["Propagator"] = {
  l1^2 - msq,
  (l1 + p1 + p2)^2 - msq,
  l2^2 - msq,
  (l2 + p1 + p2)^2 - msq,
  (l1 + p1)^2,
  (l1 - l2)^2,
  (l2 + p3)^2,
  (l2 + p1)^2,
  (l1 + p3)^2
};
AMFlowInfo["Numeric"] = {
  s -> -3,
  t -> -1,
  msq -> 1
};
AMFlowInfo["NThread"] = 4;

targets = { j[bhabha2pb, 1, 1, 1, 1, 1, 1, 1, 0, 0] };
epslist = {1/100000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "moller_scalar_qed", "mma_refs",
    "bhabha2L_PB_euclid_eps00001_mma_cache"}]];

Print["==BENCH== mma bhabha 2L planar double-box (PB), Euclidean region."];
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
