
import Base: showerror

struct InconsistentSystem <: Exception
    operator::Hecke.Generic.MatElem
    value::Hecke.Generic.MatElem
end

struct IncompatibleDimensions <: Exception    
end

function showerror(io::IO, e::InconsistentSystem)
    print(io, "Inconsistent linear system:")
    print(io, "Matrix:")
    print(io, e.operator)
    print(io, "Target:")
    print(io, e.value)
end

function showerror(io::IO, e::IncompatibleDimensions)
    print(io, "Incompatible dimensions in linear solve.")
end

