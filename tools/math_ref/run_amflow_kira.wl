(* ::Package:: *)

(* Generic AMFlow+Kira reference runner.

   Usage:
     wolfram -script run_amflow_kira.wl <config.wl>

   The config file must define the symbols documented below before
   returning. We then call AMFlow's SolveIntegrals[...] and dump
   the result as JSON into tests/data/math_ref/<name>.json.

   Required config symbols:
     name            string, output file basename (no extension)
     family          symbol, family head (e.g. bubble, box1)
     loops           list of symbols, e.g. {l1, l2}
     legs            list of symbols
     conservation    list of replacement rules, can be {}
     replacement     list of replacement rules
     propagators     list of inverse propagators
     numeric         list of replacement rules for invariants
     targets         list of integrals (j[family, ...])
     precision       integer
     epsorder        integer

   Optional:
     nthread         integer (default 4) *)

cfgFile = Environment["AMF_REF_CFG"];
If[cfgFile === $Failed || cfgFile === None,
  Print["Set environment variable AMF_REF_CFG=<config.wl>"];
  Quit[1]];
If[!FileExistsQ[cfgFile],
  Print["Config file not found: ", cfgFile];
  Quit[1]];

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
amfdir  = FileNameJoin[{current, "..", "..", "reference", "amflow-master"}];
Get[FileNameJoin[{amfdir, "AMFlow.m"}]];

(* Force Kira backend; install.m is loaded automatically, then we override paths. *)
SetReductionOptions["IBPReducer" -> "Kira"];
$KiraExecutable   = "/usr/local/bin/kira";
$FermatExecutable = "/usr/share/Ferl7/fer64";

(* User config defines the AMFlow inputs. *)
nthread = 4;
Get[cfgFile];

AMFlowInfo["Family"]       = family;
AMFlowInfo["Loop"]         = loops;
AMFlowInfo["Leg"]          = legs;
AMFlowInfo["Conservation"] = conservation;
AMFlowInfo["Replacement"]  = replacement;
AMFlowInfo["Propagator"]   = propagators;
AMFlowInfo["Numeric"]      = numeric;
AMFlowInfo["NThread"]      = nthread;
If[ValueQ[cut],          AMFlowInfo["Cut"]          = cut];
If[ValueQ[prescription], AMFlowInfo["Prescription"] = prescription];

Print["==RUNNING== ", name, " precision=", precision, " epsorder=", epsorder];

t0 = AbsoluteTime[];
sol = SolveIntegrals[targets, precision, epsorder];
t1 = AbsoluteTime[];
Print["==SOL_TIME== ", t1 - t0, "s"];
Print["==RAW=="];
Print[sol];

toLaurent[expr_] := Module[{ser, coeffs, leading, n},
  ser = Series[expr, {eps, 0, epsorder}];
  If[Head[ser] =!= SeriesData,
    Return[{{0, ToString[N[expr, precision], InputForm]}}]];
  leading = ser[[4]];
  coeffs  = ser[[3]];
  Table[{leading + n - 1, ToString[N[coeffs[[n]], precision], InputForm]},
        {n, 1, Length[coeffs]}]
];

records = Table[
  Module[{key, expr, lst},
    key  = ToString[sol[[i, 1]], InputForm];
    expr = sol[[i, 2]];
    lst  = toLaurent[expr];
    {key, lst}], {i, 1, Length[sol]}];

emitNum[s_] := "\"" <> s <> "\"";
emitTerm[t_] := "[" <> ToString[t[[1]]] <> ", " <> emitNum[t[[2]]] <> "]";
emitRecord[r_] := "    {\n      \"integral\": " <> emitNum[r[[1]]] <>
  ",\n      \"laurent\": [" <>
  StringJoin[Riffle[emitTerm /@ r[[2]], ", "]] <> "]\n    }";

emitRule[r_] := "\"" <> ToString[r[[1]], InputForm] <> "\": \"" <>
  ToString[r[[2]], InputForm] <> "\"";
numericJson = "{ " <>
  StringJoin[Riffle[emitRule /@ numeric, ", "]] <> " }";

json = "{\n  \"family\": \"" <> ToString[family] <>
  "\",\n  \"numeric\": " <> numericJson <>
  ",\n  \"precision\": " <> ToString[precision] <>
  ",\n  \"epsorder\": " <> ToString[epsorder] <>
  ",\n  \"results\": [\n" <>
  StringJoin[Riffle[emitRecord /@ records, ",\n"]] <>
  "\n  ]\n}\n";

outdir = FileNameJoin[{current, "..", "..", "tests", "data", "math_ref"}];
If[!DirectoryQ[outdir], CreateDirectory[outdir]];
outfile = FileNameJoin[{outdir, name <> ".json"}];
WriteString[outfile, json];
Print["==WROTE== ", outfile];

Quit[];
