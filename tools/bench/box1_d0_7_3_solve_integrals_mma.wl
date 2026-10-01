(* ::Package:: *)

(*
  D0 != 4 parity bench.  Same family + target as
  tools/bench/box1_single_solve_integrals_mma.wl (1-loop box at
  s = 100, t = -1, target j[box1, 1, 0, 1, 0]) but with
  GlobalOptions::D0 set to 7/3 instead of the default 4.

  Closes the natural follow-up under docs/ROADMAP.md Milestone 3:
  the C++ port now honours GlobalOptions::d0 (mirrors AMFlow.m's
  "D0" option, AMFlow.m:262), and SolveIntegrals shifts the
  internal eps grid by (4 - D0)/2 (AMFlow.m:1342, 1351); FitEps
  fits Laurent in the user-facing eps (AMFlow.m:1356).  This bench
  is the first end-to-end check that the shift drops out and the
  Laurent expansion in the user eps matches MMA at D0 != 4.

  Implementation note: SolveIntegrals takes 3 args (jints, goal,
  order); we mirror its body inline so the BlackBoxAMFlow call can
  be backed by a committed cache directory (mirrors the existing
  box1_single bench which factors SolveIntegrals the same way).
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions["IBPReducer" -> "Kira", "DeleteBlackBoxDirectory" -> False];
Kira`$KiraExecutable = "/usr/local/bin/kira";
Kira`$FermatExecutable = "/usr/share/Ferl7/fer64";

SetAMFOptions[
  "ChopPre" -> 20,
  "RationalizePre" -> 100,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5,
  "D0" -> 7/3
];

AMFlowInfo["Family"] = box1;
AMFlowInfo["Loop"] = {l};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1 - p2 - p3};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
  (p1 + p2)^2 -> s, (p1 + p3)^2 -> t
};
AMFlowInfo["Propagator"] = {
  l^2,
  (l + p1)^2,
  (l + p1 + p2)^2,
  (l + p1 + p2 + p4)^2
};
AMFlowInfo["Numeric"] = {s -> 100, t -> -1};
AMFlowInfo["NThread"] = 4;

targets = {j[box1, 1, 0, 1, 0]};
precision = 30;
epsorder = 4;

cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "box1_d0_7_3_solve_integrals_mma_cache"}]];

Print["==BENCH== mma 1L box1 SolveIntegrals at D0 = 7/3"];
Print["==CACHE== ", cache];
Print["==D0== 7/3"];
Print["==SETTINGS== precision=", precision, " eps_order=", epsorder];

total = AbsoluteTiming[
  {epslist0, workpre, xorder} = GenerateNumericalConfig[precision, epsorder];
  (* AMFlow.m:1351 — must reference AMFlow`Private`$D0, since $D0 in this
     top-level (Global`) script is a fresh undefined symbol and would leak
     into the cached epslist, causing the child kernel ODE solver to spin
     forever on a symbolic-eps system. *)
  epslist = epslist0 + (4 - AMFlow`Private`$D0)*1/2;
  Print["==CONFIG== samples=", Length[epslist], " working_pre=", workpre,
    " x_order=", xorder];
  Print["==EPSLIST0== ", epslist0 // InputForm];
  Print["==EPSLIST_INTERNAL== ", epslist // InputForm];

  Block[{$WorkingPre = workpre, $XOrder = xorder},
    blackbox = AbsoluteTiming[
      sol = BlackBoxAMFlow[targets, epslist, cache];
    ][[1]];
    (* Mirror SolveIntegrals body (AMFlow.m:1356) but resolve Loop against
       AMFlowInfo directly: AMFlow.m:192 binds Loop in the Private context,
       not at top level, so a bare ``Loop`` here would be the undefined
       Global`Loop and Length[Global`Loop] = 0 would silently flip the
       leading-order computation. *)
    leading = -2*Length[AMFlowInfo["Loop"]];
    fit = AbsoluteTiming[
      final = Thread[Keys[sol] -> N[
        Normal@Series[
          FitEps[epslist0, #, leading]& /@ Values[sol],
          {$Eps, 0, epsorder + leading}
        ],
        precision]
      ];
    ][[1]];
  ];
][[1]];

Print["==BLACKBOX_TIME== ", blackbox];
Print["==FIT_TIME== ", fit];
Print["==TOTAL_TIME== ", total];
Print["==RESULT_INPUTFORM== ", final // InputForm];
Print["==LEADING_ORDER== ", leading];

(* AMFlow.m:181 — $Eps in the AMFlow private context resolves to
   Symbol["Global`eps"], i.e. the bare lowercase ``eps`` in Global`.
   A top-level reference to ``$Eps`` here would be the unrelated
   Global`$Eps and Coefficient[expr, Global`$Eps, ord] would silently
   collapse the whole polynomial onto order=0.  Use AMFlow`Private`$Eps
   so we read the same symbol the polynomial was built on. *)
Do[
  key = Keys[final][[i]];
  expr = Values[final][[i]];
  Print["==INTEGRAL== ", key // InputForm];
  Do[
    coeff = Coefficient[expr, AMFlow`Private`$Eps, ord];
    Print["==COEFF== order=", ord,
          " :: ", N[Re[coeff], precision] // InputForm,
          " :: ", N[Im[coeff], precision] // InputForm],
    {ord, leading, epsorder + leading}
  ],
  {i, Length[Keys[final]]}
];

Quit[];
