all:
	cd libbasedir; $(MAKE) all
	cd streesrc; $(MAKE) all
	cd mm3src; $(MAKE) all

clean:
	rm -f *~
	cd libbasedir; $(MAKE) clean
	cd streesrc; $(MAKE) clean
	cd mm3src; $(MAKE) clean

mummer:
	cd libbasedir; $(MAKE) libbase.a
	cd streesrc; $(MAKE) libstree.a
	cd mm3src; $(MAKE) mummer

splintall:
	cd libbasedir; ${MAKE} splintall
	cd streesrc; ${MAKE} splintall
	cd mm3src; ${MAKE} splintall
