
.PHONY: build
build:
	cmake build
	(cd build; make)

.PHONY: format
format:
	clang-format -i src/*.cpp
	clang-format -i src/*.hpp
