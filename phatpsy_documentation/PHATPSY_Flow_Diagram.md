```mermaid
flowchart LR
phatpsy["PHATPSY\n(main.f)"]
contrl["CONTRL\n(contrl.f)"]
aloop["ALOOP\n(aloop.f)"]
bloop["BLOOP\n(bloop.f)"]
insrtd["INSRTD\n(insrtd.f)"]
atomic["ATOMIC\n(atomic.f)"]
ewmo["EWMO\n(ewmo.f)"]
excite["EXCITE\n(excite.f)"]
popout["POPOUT\n(popout.f)"]
genbc["GENBC\n(genbc.f)"]
gencgc["GENCGC\n(gencgc.f)"]
genfac["GENFAC\n(genfac.f)"]
putcgc["PUTCGC\n(putcgc.f)"]
scf["SCF\n(scf.f)"]
analys["ANALYS\n(analys.f)"]
denmat["DENMAT\n(denmat.f)"]
esca["ESCA\n(esca.f)"]
flip["FLIP\n(flip.f)"]
fokmat["FOKMAT\n(fokmat.f)"]
maxovl["MAXOVL\n(maxovl.f)"]
neword["NEWORD\n(neword.f)"]
plotv["PLOTV\n(plotv.f)"]
punchv["PUNCHV\n(punchv.f)"]
submat["ADDMAT\n(addmat.f)"]
symop["SYMOP\n(symop.f)"]
trisq["TRISQ\n(trisq.f)"]
uthu["UTHU\n(uthu.f)"]
xyzmap["XYZMAP\n(xyzmap.f)"]
symxyz["SYMXYZ\n(symxyz.f)"]
depth["DEPTH\n(depth.f)"]
oneint["ONEINT\n(oneint.f)"]
twoint["TWOINT\n(twoint.f)"]
gensab["GENSAB\n(gensab.f)"]
genvab["GENVAB\n(genvab.f)"]
jacobi["JACOBI\n(jacobi.f)"]
insert["INSERT\n(insert.f)"]
rotsab["ROTSAB\n(rotsab.f)"]
ascale["ASCALE\n(ascale.f)"]
bscale["BSCALE\n(bscale.f)"]
chkcgc["CHKCGC\n(chkcgc.f)"]
wfct["WFCT\n(wfct.f)"]
genrot["GENROT\n(genrot.f)"]
order["ORDER\n(order.f)"]
anmbnm["ANMBNM\n(anmbnm.f)"]
ndxcgc["NDXCGC\n(ndxcgc.f)"]
lowdin["LOWDIN\n(lowdin.f)"]

%% Main Flow
phatpsy --> contrl
contrl --> gencgc
contrl --> putcgc
contrl --> genbc
contrl --> genfac
contrl --> aloop
contrl --> atomic
contrl --> scf
contrl --> ewmo

%% ALOOP Flow
aloop --> bloop
aloop --> insrtd

%% BLOOP Integration - Missing Routines
bloop --> gensab
bloop --> genvab
bloop --> rotsab
bloop --> insert

%% EWMO Flow
ewmo --> excite
ewmo --> popout

%% SCF Integration
scf --> fokmat
scf --> denmat
scf --> analys
scf --> uthu
scf --> lowdin
scf --> esca
scf --> neword
scf --> plotv
scf --> punchv
scf --> submat
scf --> symop
scf --> flip
scf --> trisq
scf --> maxovl

%% ATOMIC Integration
atomic --> oneint
atomic --> twoint
atomic --> xyzmap
atomic --> symxyz

%% TWOINT Integration
twoint --> wfct

%% LOWDIN Integration
lowdin --> jacobi

%% ANMBNM Integration
gensab --> anmbnm
anmbnm --> ascale
anmbnm --> bscale

%% NDXCGD Integration
gencgc --> ndxcgc
ndxcgc --> chkcgc
ndxcgc --> order

%% SYMXYZ Integration
symxyz --> genrot

%% XYZMAP Flow
xyzmap --> depth
