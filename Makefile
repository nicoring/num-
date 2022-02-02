
.PHONY: build
build:
	cmake --build build

.PHONY: format
format:
	clang-format -i src/*.cpp
	clang-format -i src/*.hpp
	clang-format -i tests/*.cpp

.PHONY: test
test:
	(cd build; make test)
