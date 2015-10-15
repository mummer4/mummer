# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO   " $@;
.yaggo.hpp:
	$(V_YAGGO)$(YAGGO) -o $@ $<

YAGGO_BUILT = # Append all file to be built by yaggo
BUILT_SOURCES += $(YAGGO_BUILT)
noinst_HEADERS += $(YAGGO_BUILT)
DISTCLEANFILES += $(YAGGO_BUILT)
