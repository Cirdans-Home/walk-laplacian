AMGDIR=../install
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_ext -lpsb_linsolve -lamg_prec -lpsb_prec $(PSBGPULDLIBS)
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

DFSOBJS=amg_dkatz.o data_input.o
SFSOBJS=amg_df_sample.o data_input.o
CFSOBJS=
ZFSOBJS=
LINKOPT=

EXEDIR=../runs

all: amg_df_sample amg_dkatz

amg_dkatz: $(DFSOBJS)
	$(FLINK) $(LINKOPT) $(DFSOBJS) -o amg_dkatz \
	$(AMG_LIBS)  $(PSBLAS_LIBS) $(PSBGPULDLIBS) $(LDLIBS)
	/bin/mv amg_dkatz $(EXEDIR)

amg_df_sample: $(SFSOBJS)
		$(FLINK) $(LINKOPT) $(SFSOBJS) -o amg_df_sample \
		$(AMG_LIBS)  $(PSBLAS_LIBS) $(LDLIBS)
		/bin/mv amg_df_sample $(EXEDIR)

amg_df_sample.o: data_input.o

clean:
	/bin/rm -f $(DFSOBJS) $(ZFSOBJS) $(SFSOBJS)  $(CFSOBJS) \
        *$(.mod) $(EXEDIR)/amg_dkatz $(EXEDIR)/amg_df_sample
