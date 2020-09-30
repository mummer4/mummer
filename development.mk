AM_CXXFLAGS += -Werror

# Count lines of code
.PHONY: cloc
cloc:
	cloc --force-lang="Ruby,yaggo" --force-lang="make,am" --force-lang="make,mk" \
	  --exclude-dir="gtest" --not-match-d="^build-" --ignored=cloc_ignored_src_files $(srcdir) 

# Print the value of a variable
print-%:
	@echo $($*)
