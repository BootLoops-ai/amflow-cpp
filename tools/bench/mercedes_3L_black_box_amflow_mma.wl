(* ::Package:: *)

(*
  3-loop Mercedes self-energy — Phase 3 oracle (3L topology with
  internal vertex structure, not pure banana).  6 main propagators
  forming 3 outer rails (triangle A-B-C-A) and 3 inner spokes from
  central vertex O to A, B, C.  Massless internals, 2 external legs
  carrying p with p^2 = psq.  3 ISP propagators (l_i + p)^2 round
  out the 9-element SP basis (L(L+1)/2 + L*n_ext = 6 + 3 = 9).
  Numerics: psq = -3 (Euclidean), eps = 1/1000.
  Routing (Kirchhoff at vertices A, B, C, O):
    A: external p enters; rails AB, CA, spoke OA.
    B: rails AB, BC, spoke OB.
    C: external p exits; rails BC, CA, spoke OC.
    O: spokes OA + OB + OC = 0 (internal junction).
    q_AB = l1, q_BC = l2, q_CA = l3.
    q_OA = l1 - l3 - p, q_OB = l2 - l1, q_OC = l3 - l2 + p.
*)

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
  "EndingScheme" -> {"Tradition", "SingleMass"},
  "WorkingPre" -> 120, "ChopPre" -> 20, "RationalizePre" -> 100,
  "XOrder" -> 240, "ExtraXOrder" -> 280,
  "LearnXOrder" -> -1, "TestXOrder" -> 5
];

AMFlowInfo["Family"] = mercedes3L;
AMFlowInfo["Loop"] = {l1, l2, l3};
AMFlowInfo["Leg"] = {p};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p^2 -> psq};
AMFlowInfo["Propagator"] = {
  l1^2,
  l2^2,
  l3^2,
  (l1 - l3 - p)^2,
  (l2 - l1)^2,
  (l3 - l2 + p)^2,
  (l1 + p)^2,
  (l2 + p)^2,
  (l3 + p)^2
};
AMFlowInfo["Numeric"] = {psq -> -3};
AMFlowInfo["NThread"] = 4;

targets = {j[mercedes3L, 1, 1, 1, 1, 1, 1, 0, 0, 0]};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "mercedes_3L_mma_cache"}]];

Print["==BENCH== mma 3-loop Mercedes self-energy BlackBoxAMFlow"];
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
