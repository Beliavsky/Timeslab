c&INIT1
c&INIT1
        block data init1
#include 'tslabc'
	data numpr1,numpr3,numpr5/66,54,34/
	data maxnum,maxarr,maxiv,
     1       maxrv,maxchr/200000,60,60,60,60/
	data pr1lst(1) /'READ'/,
     1 pr1lst(2) /'LIST'/,
     1 pr1lst(3) /'PLOT'/,
     1 pr1lst(4) /'SAVE'/,
     1 pr1lst(5) /'PLOTSP'/,
     1 pr1lst(6) /'HELP'/,
     1 pr1lst(7) /'MACRO'/,
     1 pr1lst(8) /'PRINT'/,
     1 pr1lst(9) /'SPEAKER'/,
     1 pr1lst(10) /'EDIT'/,
     1 pr1lst(11) /'PLOTCSP'/,
     1 pr1lst(12) /'LISTSP'/,
     1 pr1lst(13) /'INFQNT'/,
     1 pr1lst(14) /'IF'/,
     1 pr1lst(15) /'HIST'/,
     1 pr1lst(16) /'CLEAN'/,
     1 pr1lst(17) /'LENGTH'/,
     1 pr1lst(18) /'PLOT2'/,
     1 pr1lst(19) /'TIME'/
	data pr1lst(20) /'DELAY'/,
     1 pr1lst(21) /'COLOR'/,
     1 pr1lst(22) /'CORR2'/,
     1 pr1lst(23) /'CROSSP'/,
     1 pr1lst(24) /'POLAR'/,
     1 pr1lst(25) /'REG'/,
     1 pr1lst(26) /'DTFORE'/,
     1 pr1lst(27) /'TEXTCOLOR'/,
     1 pr1lst(28) /'RESCREEN'/,
     1 pr1lst(29) /'MAXMIN'/,
     1 pr1lst(30) /'ARSP2'/,
     1 pr1lst(31) /'SAVESC'/,
     1 pr1lst(32) /'PRINTER'/,
     1 pr1lst(33) /'LISTM'/,
     1 pr1lst(34) /'NOTES'/,
     1 pr1lst(35) /'ARMASEL'/,
     1 pr1lst(36) /'COEFFCSD'/
	data pr1lst(37)/'FFT'/
	data pr1lst(38)/'DOUBLE'/
	data pr1lst(39)/'DTARMA'/
	data pr1lst(40)/'PRINTSEL'/
	data pr1lst(41)/'POLYROOTS'/
	data pr1lst(42)/'ROOTSPOLY'/
	data pr1lst(43)/'ARSPCB'/
	data pr1lst(44)/'ARSPPEAK'/
	data pr1lst(45)/'ARMAPRED'/
	data pr1lst(46)/'DENSITY'/
	data pr1lst(47)/'MCHOL'/
	data pr1lst(48)/'GS'/
	data pr1lst(49)/'PLOTK'/
	data pr1lst(50)/'SEASEST'/
	data pr1lst(51)/'SEASPRED'/
	data pr1lst(52)/'RECORD'/
	data pr1lst(53)/'ARCORR2'/
	data pr1lst(54)/'PLOTSIZE'/
	data pr1lst(55)/'ERASE'/
	data pr1lst(56)/'GRMENU'/
	data pr1lst(57)/'PLOT3D'/
	data pr1lst(58)/'CONTOUR'/
	data pr1lst(59)/'DEVON'/
	data pr1lst(60)/'POINTS'/
        data pr1lst(61)/'LINES'/
        data pr1lst(62)/'AXES'/
        data pr1lst(63)/'LABELS'/
        data pr1lst(64)/'XON'/
        data pr1lst(65)/'FONT'/
        data pr1lst(66)/'PANEL'/
	data pr3lst(1) /'WN'/,
     1 pr3lst(2) /'COS'/,
     1 pr3lst(3) /'SIN'/,
     1 pr3lst(4) /'DIFF'/,
     1 pr3lst(5) /'CORR'/,
     1 pr3lst(6) /'LOGE'/,
     1 pr3lst(7) /'WINDOW'/,
     1 pr3lst(8) /'SUBMNS'/,
     1 pr3lst(9) /'DTAR'/,
     1 pr3lst(10) /'ARSP'/,
     1 pr3lst(11) /'CUMSP'/,
     1 pr3lst(12) /'MASP'/,
     1 pr3lst(13) /'ARMASP'/,
     1 pr3lst(14) /'EXTRACT'/,
     1 pr3lst(15) /'ARFILT'/,
     1 pr3lst(16) /'ARDT2'/,
     1 pr3lst(17) /'EXP'/,
     1 pr3lst(18) /'BARTTEST'/,
     1 pr3lst(19) /'ARDT'/
	data pr3lst(20) /'QTEST'/,
     1 pr3lst(21) /'EXTEND'/,
     1 pr3lst(22) /'SORT'/,
     1 pr3lst(23) /'PARTAR'/,
     1 pr3lst(24) /'LINE'/,
     1 pr3lst(25) /'ARMACORR'/,
     1 pr3lst(26) /'CORRMA'/,
     1 pr3lst(27) /'MACORR'/,
     1 pr3lst(28) /'CORRARMA'/,
     1 pr3lst(29) /'MADT'/,
     1 pr3lst(30) /'ARMADT'/,
     1 pr3lst(31) /'MULTPOLY'/,
     1 pr3lst(32) /'INVPOLY'/,
     1 pr3lst(33) /'CORRAR'/,
     1 pr3lst(34) /'ARCORR'/,
     1 pr3lst(35) /'DIVSDS'/,
     1 pr3lst(36) /'ARPART'/,
     1 pr3lst(37) /'TRANS'/,
     1 pr3lst(38) /'DIST'/
	data pr3lst(39) /'POLY'/,
     1 pr3lst(40) /'DOT'/,
     1 pr3lst(41) /'MINV'/,
     1 pr3lst(42) /'MDEL'/,
     1 pr3lst(43) /'CORRAR2'/,
     1 pr3lst(44) /'MMULT'/,
     1 pr3lst(45) /'TOEPL'/
	data pr3lst(46)/'FILT'/
	data pr3lst(47)/'REVERSE'/
	data pr3lst(48)/'EIG'/
	data pr3lst(49)/'REPLACE'/
	data pr3lst(50)/'CUM'/
	data pr3lst(51)/'PARCORR'/
	data pr3lst(52)/'SWEEP'/
	data pr3lst(53)/'ABS'/
	data pr3lst(54)/'BINOM'/
	data pr5lst(1) /'QUIT'/,
     1 pr5lst(2) /'QUIT'/,
     1 pr5lst(3) /'INFO'/,
     1 pr5lst(4) /'CLEAN'/,
     1 pr5lst(5) /'LABEL'/,
     1 pr5lst(6) /'OVERON'/,
     1 pr5lst(7) /'OVEROFF'/,
     1 pr5lst(8) /'BATCHON'/,
     1 pr5lst(9) /'BATCHOFF'/,
     1 pr5lst(10) /'PSON'/,
     1 pr5lst(11) /'PSOFF'/,
     1 pr5lst(12) /'SPEAKERON'/,
     1 pr5lst(13) /'SPEAKEROFF'/,
     1 pr5lst(14) /'ECHO'/,
     1 pr5lst(15) /'DOS'/,
     1 pr5lst(16) /'PROMPTON'/,
     1 pr5lst(17) /'PROMPTOFF'/,
     1 pr5lst(18) /'PAUSE'/,
     1 pr5lst(19) /'PAGE'/
	data pr5lst(20) /'MACRO'/,
     1 pr5lst(21) /'RESTART'/,
     1 pr5lst(22) /'SINGLEON'/,
     1 pr5lst(23) /'SINGLEOFF'/,
     1 pr5lst(24) /'CLS'/,
     1 pr5lst(25) /'ABORTON'/,
     1 pr5lst(26) /'ABORTOFF'/
	data pr5lst(27)/'HELP'/
	data pr5lst(28)/'PLOTON'/
	data pr5lst(29)/'DEVOFF'/
	data pr5lst(30)/'ENDIF'/
	data pr5lst(31)/'XON'/
	data pr5lst(32)/'FIND'/
        data pr5lst(33)/'ASKON'/
        data pr5lst(34)/'ASKOFF'/
        end





