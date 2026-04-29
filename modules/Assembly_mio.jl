# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

using Memoize
using SparseArrays

"""
    initialize_assembly!(mesh::Mesh)

Initialize the assembly process for the given mesh by computing the necessary geometric quantities.

# Arguments
- `mesh::Mesh`: The mesh object for which the assembly is initialized.
"""
function initialize_assembly!(mesh::Mesh)
    get_Bk!(mesh)
    get_detBk!(mesh)
    get_invBk!(mesh)
end

########################### GLOBAL ASSEMBLER ########################### 
"""
    assemble_global(mesh::Mesh, local_assembler!)

Assemble the global stiffness matrix and force vector for the given mesh using the provided local assembler function.

# Arguments
- `mesh::Mesh`: The mesh object.
- `local_assembler!`: A function that assembles the local stiffness matrix and force vector.

# Returns
- `K::SparseMatrixCSC`: The global stiffness matrix.
- `f::Vector`: The global force vector.
"""
function assemble_global(mesh::Mesh, local_assembler!)
    ###########################################################################
    ############################ ADD CODE HERE ################################
    ########################################################################### 
end

########################################################################
########################### LOCAL ASSEMBLERS ###########################
########################################################################

########################### POISSON PROBLEM ###########################
"""
    shapef_2DLFE(quadrule::TriQuad)

Compute the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `shapef`: The shape functions evaluated at the quadrature points.
"""
#@memoize function shapef_2DLFE(quadrule::TriQuad)
function shapef_2DLFE(quadrule::TriQuad)
    ###########################################################################
    points = quadrule.points
    x = points[1,:]
    y = points[2,:]

    phi = [
        (1 .- x .- y),          # phi_1
        x,                      # phi_2
        y                       # phi_3
    ]

    return phi

    ########################################################################### 
end

"""
    ∇shapef_2DLFE(quadrule::TriQuad)

Compute the gradients of the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `∇shapef`: The gradients of the shape functions evaluated at the quadrature points.
"""
#@memoize function ∇shapef_2DLFE(quadrule::TriQuad)
function ∇shapef_2DLFE(quadrule::TriQuad)
    ###########################################################################
    mat = [-1  1  0;
            1  0  1]   # matrice 2×3

    q = size(quadrule.points, 2)

    return repeat(mat, 1, 1, q)
    ########################################################################### 
end

"""
    poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)

Assemble the local stiffness matrix and force vector for the Poisson problem.

# Arguments
- `Ke::Matrix`: The local stiffness matrix to be assembled.
- `fe::Vector`: The local force vector to be assembled.
- `mesh::Mesh`: The mesh object.
- `cell_index::Integer`: The index of the current cell.
- `f`: The source term function.

# Returns
- `Ke`: The assembled local stiffness matrix.
- `fe`: The assembled local force vector.
"""
function poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)
    ###########################################################################
    ke = zeros(3,3)
    fe = zeros(3)
    quadrule = Q2_ref
    ak = mesh.ak[:, cell_index]
    Bk = mesh.Bk[:, :, cell_index]
    invBk = mesh.invBk[:, :, cell_index]
    pe = Bk * quadrule.points + ak
    phie = shapef_2DLFE(quadrule)
    gradphie = mult(invBk.T,∇shapef_2DLF(quadrule), axis = 1)
    for p = 1:q
        wp = quadrule[p].abs(mesh.detBk[cell_index])
        for i = 1:3
            v = phie[i,p]
            gradv = gradphie[:, i, p]
            fe[i] += wp*v*f(pe[:, p])
            for j = 1:3
                gradu = gradphie[:, j, p]
                ke[j, i] += wp * gradv *gradu
            end
        end
    end
    ########################################################################### 
end