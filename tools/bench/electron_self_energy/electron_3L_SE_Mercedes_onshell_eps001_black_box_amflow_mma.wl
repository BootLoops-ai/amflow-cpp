(* ::Package:: *)

(*
  3-loop QED electron self-energy — Mercedes-Benz rim 2-leg topology.
  Graph: central hub O + 3 photon spokes (O-A, O-B, O-C) + 5 electron
  rim arcs (A-M, M-B, B-C, C-N, N-A) where M, N are the 2 external
  attachment points. ext_in -> M, ext_out -> N, both carrying p1.
  Mass: 3 spokes = photon (massless), 5 rim arcs = electron (mass m_e).
  Loop momenta: l1 = e_OA, l2 = e_OB, l3 = e_NA.
  ON-SHELL kinematics: p1^2 = ssq = 1 = msq (electron on mass shell).
  eps = 1/1000.
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

AMFlowInfo["Family"] = eseMER;
AMFlowInfo["Loop"] = {l1, l2, l3};
AMFlowInfo["Leg"] = {p1, p2};
AMFlowInfo["Conservation"] = {p2 -> -p1};
AMFlowInfo["Replacement"] = {p1^2 -> ssq};
AMFlowInfo["Propagator"] = {
  l1^2,
  l2^2,
  (l1 + l2)^2,
  (l1 + l3)^2 - msq,
  (p1 + l1 + l3)^2 - msq,
  (p1 + l1 + l2 + l3)^2 - msq,
  (p1 + l3)^2 - msq,
  l3^2 - msq,
  (p1 + l2)^2
};
AMFlowInfo["Numeric"] = {ssq -> 1, msq -> 1};
AMFlowInfo["NThread"] = 4;

targets = { j[eseMER, 1, 1, 1, 1, 1, 1, 1, 1, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "electron_self_energy", "mma_refs",
    "electron_3L_SE_Mercedes_onshell_eps001_mma_cache"}]];

Print["==BENCH== mma electron 3L SE Mercedes rim 2-leg ON-SHELL (ssq=msq=1)."];
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
