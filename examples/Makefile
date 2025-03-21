all:
	mkdir output
	@$(MAKE) get_cms_od -C ..
	@printf "\n"
	@$(MAKE) penc_example -C ..
	@printf "\n"
	@$(MAKE) re3c_example -C ..
	@printf "\n"
	@chmod +x visualize_pencs.py
	@chmod +x visualize_re3c.py

penc_example:
	./write_penc --weights 1 4 9 49 99
	./visualize_pencs.py
	open output/penc_example.png

re3c_example:
	./write_re3c
	# ./visualize_re3cs.py
	# open output/re3c_example_.png
