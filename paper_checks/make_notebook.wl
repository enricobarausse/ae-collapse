(* Build check_equations.nb from check_equations.wl (same directory): a Title cell, the script's header comment as a
   Text cell, a Text cell on how to evaluate, then for each "(* ::Section:: *)" of the script a Section cell, a Text
   cell with the section's leading comment, and ONE Input cell with the section's code.
   Run:  wolframscript -file make_notebook.wl
   Regenerate and commit check_equations.nb after every edit of check_equations.wl (the .nb is not rebuilt automatically). *)
dir = DirectoryName[$InputFileName];
src = Import[FileNameJoin[{dir, "check_equations.wl"}], "Text"];
parts = StringSplit[src, "(* ::Section:: *)"];
secs = Rest[parts];
titleOf[s_] := First[StringCases[s, StartOfString ~~ Whitespace ... ~~ "(*" ~~ t : Shortest[__] ~~ "*)" :> t], "?"];
bodyOf[s_] := StringTrim[StringReplace[s, StartOfString ~~ Whitespace ... ~~ "(*" ~~ Shortest[__] ~~ "*)" -> "", 1]];
textOf[b_] := First[StringCases[b, StartOfString ~~ "(*" ~~ t : Shortest[__] ~~ "*)" :> StringTrim[t]], ""];
(* the header: the "(* ::Text:: *)" block of the script, one "(*...*)" line per paragraph line *)
headerText = StringRiffle[StringCases[Last[StringSplit[First[parts], "(* ::Text:: *)"]],
    StartOfLine ~~ "(*" ~~ t : Shortest[___] ~~ "*)" ~~ EndOfLine :> t], "\n"];
cells = {Cell["Equation checks for \"Gravitational collapse in Lorentz-violating gravity and the spin-0 mode\"", "Title"],
   Cell[headerText, "Text"],
   Cell["This notebook is generated from check_equations.wl by make_notebook.wl (same code, one Input cell per section). Evaluate it top to bottom (Evaluation > Evaluate Notebook); about 40 s with Mathematica 14.2. Each section prints one PASS/FAIL/CAVEAT/SKIP line; the last cell prints the summary (with the list of CAVEATs). For the snippet cross-check S81 set the environment variable PAPER_TEX to the path of main.tex, or put main.tex in the kernel's current directory; otherwise S81 reports SKIP. See README_mathematica.md for a run log and the equation list.", "Text"]};
Do[With[{tt = titleOf[s], b = bodyOf[s]},
   cells = Join[cells, {Cell[tt, "Section"]},
     If[textOf[b] =!= "", {Cell[textOf[b], "Text"]}, {}],
     {Cell[b, "Input"]}]], {s, secs}];
Export[FileNameJoin[{dir, "check_equations.nb"}], Notebook[cells]];
Print["check_equations.nb: ", Length[cells], " cells, ", Length[secs], " Input cells"];
