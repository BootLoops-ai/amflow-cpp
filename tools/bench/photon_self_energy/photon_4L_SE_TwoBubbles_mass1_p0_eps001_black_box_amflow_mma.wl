(* ::Package:: *)

(*
  4-loop photon self-energy TwoBubbles topology — staged mass1.
  Only propagator D1 = l1^2 carries mass me^2; all other 7 originally
  electron edges (D2, D3, D5, D6, D7, D8, D11) are made massless. The
  3 photon edges (D4, D9, D10) remain massless. ISPs are D12-D14.
  Kinematics: p1^2 = ssq = 0 (massless ext), me^2 = msq = 1, eps = 1/1000.
  ibp_dot = 5.
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", "..", ".."}];
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
  "EndingScheme" -> {"Tradition", "SingleMass"},
  "WorkingPre" -> 120, "ChopPre" -> 20, "RationalizePre" -> 100,
  "XOrder" -> 240, "ExtraXOrder" -> 280,
  "LearnXOrder" -> -1, "TestXOrder" -> 5
];

AMFlowInfo["Family"] = pseTB1;
AMFlowInfo["Loop"] = {l1, l2, l3, l4};
AMFlowInfo["Leg"] = {p1, p2};
AMFlowInfo["Conservation"] = {p2 -> -p1};
AMFlowInfo["Replacement"] = {p1^2 -> ssq};
AMFlowInfo["Propagator"] = {
  l1^2 - msq,
  (l1 + p1)^2,
  l2^2,
  (l1 + l2 + p1)^2,
  l3^2,
  (l3 + p1)^2,
  (l1 + l2 + l3 + p1)^2,
  l4^2,
  (l3 + l4 + p1)^2,
  (l1 + l2 + l3 + l4 + p1)^2,
  (l1 + l3 + l4 + p1)^2,
  (p1 + l2)^2,
  (p1 + l4)^2,
  (l1 + l3)^2
};
AMFlowInfo["Numeric"] = {ssq -> 0, msq -> 1};
AMFlowInfo["NThread"] = 4;

targets = { j[pseTB1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "photon_self_energy", "mma_refs",
    "photon_4L_SE_TwoBubbles_mass1_p0_eps001_mma_cache"}]];

Print["==BENCH== mma photon 4L SE TwoBubbles mass1 (p^2=0, msq=1, dot=5)."];
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
