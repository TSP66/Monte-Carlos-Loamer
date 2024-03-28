using MethodAnalysis
using JuliaInterpreter

visit(Base) do item
    isa(item, Module) && push!(JuliaInterpreter.compiled_modules, item)
        true
end

print("Debugging")