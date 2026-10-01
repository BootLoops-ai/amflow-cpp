(* ::Package:: *)

(*
  4-loop massless cutbanana — Phase 3 oracle (multi-cut Cutkosky axis).
  All 5 internal lines of the 4-loop equal-mass banana topology cut
  simultaneously, corresponding to a 5-particle Cutkosky cut.  Extends
  the cut series (cutbubble: 2-cut, cutsunrise: 3-cut, cutbanana_3L:
  4-cut) to 5 simultaneous on-shell propagators.  Target is the corner
  j[cutbanana4, 1,1,1,1,1, 0,...,0] at s = 4, eps = 1/100.
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

(* BlackBoxDot = 5 because the top sector indices are all 1 (JDot = 0);
   AMFlow.m's Max[$BlackBoxDot, JDot/@...] cannot raise it on its own.
   Same constraint as cutbanana_3L and banana_3loop. *)
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
  "EndingScheme" -> {"Tradition", "Cutkosky", "SingleMass"},
  "WorkingPre" -> 120,
  "ChopPre" -> 20,
  "RationalizePre" -> 100,
  "XOrder" -> 200,
  "ExtraXOrder" -> 240,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = cutbanana4;
AMFlowInfo["Loop"] = {l1, l2, l3, l4};
AMFlowInfo["Leg"] = {p};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p^2 -> s};
AMFlowInfo["Propagator"] = {
  l1^2,
  l2^2,
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
AMFlowInfo["Cut"] = {1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
AMFlowInfo["Numeric"] = {s -> 4};
AMFlowInfo["NThread"] = 4;

targets = {
  j[cutbanana4, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
};
epslist = {1/100};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "cutbanana_4L_eps001_mma_cache"}]];

Print["==BENCH== mma 4-loop massless cutbanana BlackBoxAMFlow sampled"];
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
