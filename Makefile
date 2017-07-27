.PHONY: clean
clean:
	rm *.itp* *.gro* *.dat* *.json* *.xtc* | true
	rm -r *.ff | true

.PHONY: test
test:
	py.test test/

.PHONY: check
check: test
