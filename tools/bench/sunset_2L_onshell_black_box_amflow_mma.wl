(* ::Package:: *)

(*
  2-loop equal-mass sunrise (sunset) with EXTERNAL ON-SHELL: s = msq,
  i.e. p^2 = m^2.  Sits between Euclidean (s < 0, banana_3loop) and
  threshold (s = 9m^2, sunrise_2L_threshold).  Distinctive: on-shell
  external kinematics — a common case in self-energy diagrams.
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions[
  "IBPReducer" -> "Kira",
  "BlackBoxRank" -> 0,
  "BlackBoxDot" -> 1,
  "DeleteBlackBoxDirectory" -> False
];
Kira`$KiraExecutable = "/usr/local/bin/kira";
Kira`$FermatExecutable = "/usr/share/Ferl7/fer64";

SetAMFOptions[
  "AMFMode" -> {"Prescription", "Mass", "Propagator"},
  "EndingScheme" -> {"Tradition", "SingleMass"},
  "WorkingPre" -> 120,
  "ChopPre" -> 20,
  "RationalizePre" -> 100,
  "XOrder" -> 200,
  "ExtraXOrder" -> 240,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = sunsetos;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p^2 -> s};
AMFlowInfo["Propagator"] = {
  l1^2 - msq,
  l2^2 - msq,
  (l1 + l2 - p)^2 - msq,
  (l1 + p)^2,
  (l2 + p)^2
};
AMFlowInfo["Numeric"] = {s -> 1, msq -> 1};
AMFlowInfo["NThread"] = 4;

targets = {j[sunsetos, 1, 1, 1, 0, 0]};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "sunset_2L_onshell_mma_cache"}]];

Print["==BENCH== mma 2-loop sunset on-shell BlackBoxAMFlow"];
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
    N[Re[val], 40] // InputForm, " :: ", N[Im[val], 40] // InputForm],
  {i, Length[Keys[sol]]}
];

Quit[];
