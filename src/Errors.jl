
import Base: showerror

struct InconsistentSystemError <: Exception
    operator::Hecke.Generic.MatElem
    value::Hecke.Generic.MatElem
end

struct IncompatibleOptionsError <: Exception
    msg
end

struct ConvergenceFailureError <: Exception
    msg
end

struct InsufficientPrecisionError <: Exception end

##############################################################################################
#
#   Show error functions
#
##############################################################################################

function showerror(io::IO, e::InconsistentSystemError)
    print(io, "Inconsistent linear system:\n")
    print(io, "Matrix:\n")
    print(io, e.operator)
    print(io, "\n\nTarget:\n")
    print(io, e.value)
end

function showerror(io::IO, e::IncompatibleOptionsError)
    print(io, "IncompatibleOptions: ", e.msg)
end

function showerror(io::IO, e::ConvergenceFailureError)
    print(io, "Convergence Failure: ", e.msg)
end
