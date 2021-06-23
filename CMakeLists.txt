cmake_minimum_required(VERSION 3.10)
project(mdf)

## Essa flag deixa o codigo mais rapido
set(CMAKE_CXX_FLAGS "-O2")

## Voce pode usar o comando file para obter a lista de todos
## arquivos na raiz do projeto .cpp e de dentro do diretorio src
file(GLOB source_files *.cpp src/*.cpp)

# Com esse comando, vc consegue indicar para o compilador
# a localização dos seus arquivos .h
include_directories(include)

# Observe abaixo o uso da variável source_files
add_executable(mdf.exe ${source_files})