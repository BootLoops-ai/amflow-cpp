(* ::Package:: *)

(*
  3-loop QED electron self-energy — Transverse-Cross (TX).
  Pairing (1,6)(2,4)(3,5): outer envelope photon gamma_16, plus inner
  crossed photons gamma_24 and gamma_35.

  TOPOLOGICAL NOTE: TX is 2-particle reducible (cutting e2 + e6 splits
  the graph). e2 (v1→v2) and e6 (v5→v6) carry the SAME momentum (p - l1)
  in any loop routing — geometric fact, not artifact. Hence the IBP
  family has only 7 distinct main propagators + 2 ISPs (instead of 8+1),
  and the target integral has index 2 on the shared propagator
  position #1: [2, 1, 1, 1, 1, 1, 1, 0, 0]. This is a valid test of
  AMFlow with a doubled propagator / shared-momentum graph.

  Kinematics: p1^2 = ssq = -3, m^2 = 1, eps = 1/1000.
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
  "WorkingPre" -> 120, "ChopPre" -> 20, "RationalizePre" -> 100,
  "XOrder" -> 240, "ExtraXOrder" -> 280,
  "LearnXOrder" -> -1, "TestXOrder" -> 5
];

AMFlowInfo["Family"] = eseTX;
AMFlowInfo["Loop"] = {l1, l2, l3};
AMFlowInfo["Leg"] = {p1, p2};
AMFlowInfo["Conservation"] = {p2 -> -p1};
AMFlowInfo["Replacement"] = {p1^2 -> ssq};
AMFlowInfo["Propagator"] = {
  (p1 - l1)^2 - msq,
  (p1 - l1 - l2)^2 - msq,
  (p1 - l1 - l2 - l3)^2 - msq,
  (p1 - l1 - l3)^2 - msq,
  l1^2,
  l2^2,
  l3^2,
  (p1 - l2)^2,
  (p1 - l3)^2
};
AMFlowInfo["Numeric"] = {ssq -> -3, msq -> 1};
AMFlowInfo["NThread"] = 4;

targets = { j[eseTX, 2, 1, 1, 1, 1, 1, 1, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "electron_self_energy", "mma_refs",
    "electron_3L_SE_TX_euclid_eps001_mma_cache"}]];

Print["==BENCH== mma electron 3L SE TX (1,6)(2,4)(3,5) euclid (ssq=-3, msq=1)."];
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
