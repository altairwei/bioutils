from conans import ConanFile, CMake, tools

class BioUtilsConan(ConanFile):
    name = "BioUtils"
    version = "0.0.1"
    license = "MIT"
    url = "https://github.com/altairwei/bioutils.git"
    description = "Bioutils includes most of the basic command-line tools that are expected in bioinformatics."
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"
    requires = "glib/2.58.3@bincrafters/stable"
    
    def build(self):
        cmake = CMake(self)
        cmake.definitions["CONAN_INSTALL_MANUALLY"] = "ON"
        cmake.configure()
        cmake.build()