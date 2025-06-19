
mutable struct Control
    u::Float64
end

# Getter
get_u(c::Control) = c.u

# Setter
function set_u!(c::Control, u::Float64)
    c.u = u
    return c
end