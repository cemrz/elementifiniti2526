# Author: Ivan Bioli (https://github.com/IvanBioli)

"""
    triarea(V1, V2, V3)

Calculate the area of a triangle given its vertices.

# Arguments
- `V1`: The first vertex of the triangle.
- `V2`: The second vertex of the triangle.
- `V3`: The third vertex of the triangle.

# Returns
- `area::Float64`: The area of the triangle.
"""
function triarea(V1, V2, V3)
    ###########################################################################
    return 1/2 * det([V1 V2 V3; ones(1,3)])
    ###########################################################################
end

"""
    Q0(p, T, u)

Perform numerical integration using the Q0 quadrature rule (i.e., baricenter formula) over a mesh.
This quadrature rule has order 1.

# Arguments
- `p::Matrix`: The coordinates of the mesh nodes.
- `T::Matrix`: The connectivity matrix of the mesh elements.
- `u::Function`: The function to be integrated.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""
function Q0(p, T, u)
    ###########################################################################
    Q0 = 0.0
    for j = 1:size(T, 2)
        v1 = p[:,T[1,j]]
        v2 = p[:,T[2,j]]
        v3 = p[:,T[3,j]]
        b_T= (v1 + v2 + v3) ./ 3
        area = triarea(v1, v2, v3)
        Q0 += area * u(b_T)
    end
    ###########################################################################
    return Q0
end

"""
    Q1(p, T, u)

Perform numerical integration using the Q1 quadrature rule (i.e., vertex formula) over a mesh.
This quadrature rule has order 1.

# Arguments
- `p::Matrix`: The coordinates of the mesh nodes.
- `T::Matrix`: The connectivity matrix of the mesh elements.
- `u::Function`: The function to be integrated.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""
function Q1(p, T, u)
    ###########################################################################
    Q1 = 0.0
    for j = 1:size(T, 2)
        v1 = p[:,T[1,j]]
        v2 = p[:,T[2,j]]
        v3 = p[:,T[3,j]]
        temp = (u(v1) + u(v2) + u(v3)) ./ 3
        area = triarea(v1, v2, v3)
        Q1 += area * temp
    end
    ###########################################################################
    return Q1
end

"""
    Q2(p, T, u)

Perform numerical integration using the Q2 quadrature rule (i.e., midpoints rule) over a mesh.
This quadrature rule has order 2.

# Arguments
- `p::Matrix`: The coordinates of the mesh nodes.
- `T::Matrix`: The connectivity matrix of the mesh elements.
- `u::Function`: The function to be integrated.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""
function Q2(p, T, u)
    ###########################################################################
    Q2 = 0.0
    for j = 1:size(T, 2)
        v1 = p[:,T[1,j]]
        v2 = p[:,T[2,j]]
        v3 = p[:,T[3,j]]
        m1 = (v2 + v3) ./ 2
        m2 = (v1 + v3) ./ 2 
        m3 = (v1 + v2) ./ 2
        temp = (u(m1) + u(m2) + u(m3)) ./ 3
        area = triarea(v1, v2, v3)
        Q2 += area * temp
        end
    ###########################################################################
    return Q2
end