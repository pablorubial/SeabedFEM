module MarineForm

# Loading the packages
using Distributions
using Gmsh

export generate_animals, PorousCircle, RigidCircle, PorousClamp, RigidClamp, RigidEllipse, generate_geometry

"Definition of the abstract type MarineForms."
abstract type MarineForms end

"Definition of the abstract type Circle."
abstract type Circle <: MarineForms end

"Definition of the abstract type Clamp."
abstract type Clamp <: MarineForms end

"Definition of the abstract type Ellipse"
abstract type Ellipse <: MarineForms end


"Definition of the Porous Circle object."
struct PorousCircle <: Circle
    x::Float64 # x-coordinate of the center of the circle
    y::Float64 # y-coordinate of the center of the circle
    re::Float64 # Radius of the circle
    tol::Float64 # Tolerance of the circle to avoid collisions with the boundaries of the physical domain
end

"Definition of the Rigid Circle object."
struct RigidCircle <: Circle
    x::Float64 # x-coordinate of the center of the circle
    y::Float64 # y-coordinate of the center of the circle
    re::Float64 # Radius of the circle
    tol::Float64 # Tolerance of the circle to avoid collisions with the boundaries of the physical domain
end

"Definition of the Porous Clamp object."
struct PorousClamp <: Clamp
    x::Float64 # x-coordinate of the center of the clamp
    y::Float64 # y-coordinate of the center of the clamp
    θo::Float64 # Angle of the opening of the clamp
    θb::Float64 # Angle of the bisector of the opening of the clamp
    ri::Float64 # Inner radius of the clamp
    re::Float64 # Outer radius of the clamp 
    tol::Float64 # Tolerance of the clamp to avoid collisions with the boundaries of the physical domain
end

"Definition of the Rigid Clamp object."
struct RigidClamp <: Clamp
    x::Float64 # x-coordinate of the center of the clamp
    y::Float64 # y-coordinate of the center of the clamp
    θo::Float64 # Angle of the opening of the clamp
    θb::Float64 # Angle of the bisector of the opening of the clamp
    ri::Float64 # Inner radius of the clamp
    re::Float64 # Outer radius of the clamp
    tol::Float64 # Tolerance of the clamp
end

"Definition of the Rigid Elliptical Clamp object."

struct RigidEllipse <: Ellipse
    x::Float64 # x-coordinate of the center of the ellipse
    y::Float64 # y-coordinate of the center of the ellipse
    a::Float64 # Semi-major axis
    b::Float64 # Semi-minor axis
    e::Float64 # Thickness of the elliptical clamp
    θ₀::Float64 # Angle of the opening of the elliptical clamp
    α::Float64 # Angle of orientation of the elliptical clamp
    tol::Float64 # Tolerance of the elliptical clamp
end

x(r, θ) = r * cos(θ) # x-coordinate of the point (r, θ)
y(r, θ) = r * sin(θ) # y-coordinate of the point (r, θ)
R(α) = [cos(α) -sin(α); sin(α) cos(α)] # Rigid-body rotation matrix


"This function generates a Circle object based on the given parameters.
The first argument is the type of the object to generate, and the second argument
is a dictionary with the parameters of the object."


function generate(::Type{RigidCircle}, params::Dict{Symbol,Any})
    x_range = get(params, :x_range, (0.0, 1.0))
    y_range = get(params, :y_range, (0.0, 1.0))
    r_distribution = get(params, :r_distribution, Normal(1.0, 0.2))
    
    # Needed to compute the tolerance
    r = rand(r_distribution)
    
    return RigidCircle(
        x_range[1] + (x_range[2] - x_range[1]) * rand(), 
        y_range[1] + (y_range[2] - y_range[1]) * rand(),
        r,
        0.5 * r
    )
end

function generate(::Type{PorousCircle}, params::Dict{Symbol,Any})
    x_range = get(params, :x_range, (0.0, 1.0))
    y_range = get(params, :y_range, (0.0, 1.0))
    r_distribution = get(params, :r_distribution, Normal(1.0, 0.2))
    
    # Needed to compute the tolerance
    r = rand(r_distribution)
    
    return PorousCircle(
        x_range[1] + (x_range[2] - x_range[1]) * rand(), 
        y_range[1] + (y_range[2] - y_range[1]) * rand(),
        r,
        0.5 * r
    )
end

"This function generates a Clamp object based on the given parameters.
The first argument is the type of the object to generate, and the second argument
is a dictionary with the parameters of the object.
"
function generate(::Type{RigidClamp}, params::Dict{Symbol,Any})
    x_range = get(params, :x_range, (0.0, 1.0))
    y_range = get(params, :y_range, (0.0, 1.0))
    θo_range = get(params, :θo_range, (0.0, 2π))
    θb_range = get(params, :θb_range, (0.0, 2π))
    ri_distribution = get(params, :ri_distribution, Normal(0.5, 0.1))
    re_distribution = get(params, :re_distribution, Normal(1.5, 0.2))
    
    # Needed to compute the tolerance
    ri = rand(ri_distribution)
    re = rand(re_distribution)

    return RigidClamp(
        x_range[1] + (x_range[2] - x_range[1]) * rand(), 
        y_range[1] + (y_range[2] - y_range[1]) * rand(),
        θo_range[1] + (θo_range[2] - θo_range[1]) * rand(),
        θb_range[1] + (θb_range[2] - θb_range[1]) * rand(), 
        ri,
        re,
        0.5 * (re - ri)
    )
end

function generate(::Type{PorousClamp}, params::Dict{Symbol,Any})
    x_range = get(params, :x_range, (0.0, 1.0))
    y_range = get(params, :y_range, (0.0, 1.0))
    θo_range = get(params, :θo_range, (0.0, 2π))
    θb_range = get(params, :θb_range, (0.0, 2π))
    ri_distribution = get(params, :ri_distribution, Normal(0.5, 0.1))
    re_distribution = get(params, :re_distribution, Normal(1.5, 0.2))
    
    # Needed to compute the tolerance
    ri = rand(ri_distribution)
    re = rand(re_distribution)
    
    return PorousClamp(
        x_range[1] + (x_range[2] - x_range[1]) * rand(), 
        y_range[1] + (y_range[2] - y_range[1]) * rand(),
        θo_range[1] + (θo_range[2] - θo_range[1]) * rand(),
        θb_range[1] + (θb_range[2] - θb_range[1]) * rand(), 
        ri,
        re,
        0.5 * (re - ri)
    )
end

function generate(::Type{RigidEllipse}, params::Dict{Symbol,Any})
    x_range = get(params, :x_range, (0.0, 1.0))
    y_range = get(params, :y_range, (0.0, 1.0))
    a_distribution = get(params, :a_distribution, Normal(1.0, 0.2))
    b_distribution = get(params, :b_distribution, Normal(0.5, 0.1))
    e_range = get(params, :e_range, (0.01, 0.02))
    θo_range = get(params, :θo_range, (0.0, 2π))
    α_range = get(params, :α_range, (0.0, 2π))
    
    # Needed to compute the tolerance
    a = rand(a_distribution)
    b = rand(b_distribution)
    e = rand(e_range)
    
    return RigidEllipse(
        x_range[1] + (x_range[2] - x_range[1]) * rand(), 
        y_range[1] + (y_range[2] - y_range[1]) * rand(),
        a,
        b,
        e_range[1] + (e_range[2] - e_range[1]) * rand(),
        θo_range[1] + (θo_range[2] - θo_range[1]) * rand(),
        α_range[1] + (α_range[2] - α_range[1]) * rand(),
        0.5 * (a + e)
    )
    
end

function generate_geometry(animal::PorousCircle, h::Float64)
    
    # Create points that conforms the circle
    p1 = gmsh.model.geo.addPoint(animal.x + animal.re, animal.y, 0, h)
    p2 = gmsh.model.geo.addPoint(animal.x, animal.y, 0, h)
    p3 = gmsh.model.geo.addPoint(animal.x - animal.re, animal.y, 0, h)

    # Create arcs
    l1 = gmsh.model.geo.addCircleArc(p1, p2, p3) # Rock_circle
    l2 = gmsh.model.geo.addCircleArc(p3, p2, p1) # Rock_circle
    
    # Create the curve loop
    cl = gmsh.model.geo.addCurveLoop([l1, l2])

    return cl, [l1, l2]

end

function generate_geometry(animal::RigidCircle, h::Float64)
    
    # Create points that conforms the circle
    p1 = gmsh.model.geo.addPoint(animal.x + animal.re, animal.y, 0, h)
    p2 = gmsh.model.geo.addPoint(animal.x, animal.y, 0, h)
    p3 = gmsh.model.geo.addPoint(animal.x - animal.re, animal.y, 0, h)

    # Create arcs
    l1 = gmsh.model.geo.addCircleArc(p1, p2, p3) # Rock_circle
    l2 = gmsh.model.geo.addCircleArc(p3, p2, p1) # Rock_circle
    
    # Create the curve loop
    cl = gmsh.model.geo.addCurveLoop([l1, l2])

    return cl, [l1, l2]
    
end

function generate_geometry(animal::RigidClamp, h::Float64)
    
    centerPoint = gmsh.model.geo.addPoint(animal.x, animal.y, 0, h)

    p1= gmsh.model.geo.addPoint(animal.x + x(animal.ri, animal.θb + animal.θo/2), animal.y + y(animal.ri, animal.θb + animal.θo/2), 0, h)
    p2 = gmsh.model.geo.addPoint(animal.x + x(animal.ri, animal.θb + animal.θo/2 + pi/2), animal.y + y(animal.ri, animal.θb + animal.θo/2 + pi/2), 0, h)
    p3 = gmsh.model.geo.addPoint(animal.x + x(animal.ri, animal.θb + animal.θo/2 + pi), animal.y + y(animal.ri, animal.θb + animal.θo/2 + pi), 0, h)
    p4 = gmsh.model.geo.addPoint(animal.x + x(animal.ri, animal.θb - animal.θo/2), animal.y + y(animal.ri, animal.θb - animal.θo/2), 0, h)

    p5 = gmsh.model.geo.addPoint(animal.x + x(animal.re, animal.θb + animal.θo/2), animal.y + y(animal.re, animal.θb + animal.θo/2), 0, h)
    p6 = gmsh.model.geo.addPoint(animal.x + x(animal.re, animal.θb + animal.θo/2 + pi/2), animal.y + y(animal.re, animal.θb + animal.θo/2 + pi/2), 0, h)
    p7 = gmsh.model.geo.addPoint(animal.x + x(animal.re, animal.θb + animal.θo/2 + pi), animal.y + y(animal.re, animal.θb + animal.θo/2 + pi), 0, h)
    p8 = gmsh.model.geo.addPoint(animal.x + x(animal.re, animal.θb - animal.θo/2), animal.y + y(animal.re, animal.θb - animal.θo/2), 0, h)

    # Create lines and arcs
    arc_i1 = gmsh.model.geo.addCircleArc(p1, centerPoint, p2) # Curve 5
    arc_i2 = gmsh.model.geo.addCircleArc(p2, centerPoint, p3)  # curve 6
    arc_i3 = gmsh.model.geo.addCircleArc(p3, centerPoint, p4)  # curve 6

    line_1 = gmsh.model.geo.addLine(p1, p5) # curve 7

    arc_e1 = gmsh.model.geo.addCircleArc(p5, centerPoint, p6) # curve 8
    arc_e2 = gmsh.model.geo.addCircleArc(p6, centerPoint, p7) # curve 9
    arc_e3 = gmsh.model.geo.addCircleArc(p7, centerPoint, p8) # curve 10

    line_2 = gmsh.model.geo.addLine(p8, p4) # Curve 10

    # Create curve loops with the lines and arcs
    cl = gmsh.model.geo.addCurveLoop([-arc_i3, -arc_i2, -arc_i1, line_1, arc_e1, arc_e2, arc_e3, line_2])

    return cl, [arc_i1, arc_i2, arc_i3, line_1, arc_e1, arc_e2, arc_e3, line_2]  
end

function generate_geometry(animal::PorousClamp, h::Float64)
    centerPoint = gmsh.model.geo.addPoint(animal.x, animal.y, 0, h)

    p1= gmsh.model.geo.addPoint(animal.x + x(animal.ri, animal.θb + animal.θo/2), animal.y + y(animal.ri, animal.θb + animal.θo/2), 0, h)
    p2 = gmsh.model.geo.addPoint(animal.x + x(animal.ri, animal.θb + animal.θo/2 + pi/2), animal.y + y(animal.ri, animal.θb + animal.θo/2 + pi/2), 0, h)
    p3 = gmsh.model.geo.addPoint(animal.x + x(animal.ri, animal.θb + animal.θo/2 + pi), animal.y + y(animal.ri, animal.θb + animal.θo/2 + pi), 0, h)
    p4 = gmsh.model.geo.addPoint(animal.x + x(animal.ri, animal.θb - animal.θo/2), animal.y + y(animal.ri, animal.θb - animal.θo/2), 0, h)

    p5 = gmsh.model.geo.addPoint(animal.x + x(animal.re, animal.θb + animal.θo/2), animal.y + y(animal.re, animal.θb + animal.θo/2), 0, h)
    p6 = gmsh.model.geo.addPoint(animal.x + x(animal.re, animal.θb + animal.θo/2 + pi/2), animal.y + y(animal.re, animal.θb + animal.θo/2 + pi/2), 0, h)
    p7 = gmsh.model.geo.addPoint(animal.x + x(animal.re, animal.θb + animal.θo/2 + pi), animal.y + y(animal.re, animal.θb + animal.θo/2 + pi), 0, h)
    p8 = gmsh.model.geo.addPoint(animal.x + x(animal.re, animal.θb - animal.θo/2), animal.y + y(animal.re, animal.θb - animal.θo/2), 0, h)

    # Create lines and arcs
    arc_i1 = gmsh.model.geo.addCircleArc(p1, centerPoint, p2) # Curve 5
    arc_i2 = gmsh.model.geo.addCircleArc(p2, centerPoint, p3)  # curve 6
    arc_i3 = gmsh.model.geo.addCircleArc(p3, centerPoint, p4)  # curve 6

    line_1 = gmsh.model.geo.addLine(p1, p5) # curve 7

    arc_e1 = gmsh.model.geo.addCircleArc(p5, centerPoint, p6) # curve 8
    arc_e2 = gmsh.model.geo.addCircleArc(p6, centerPoint, p7) # curve 9
    arc_e3 = gmsh.model.geo.addCircleArc(p7, centerPoint, p8) # curve 10

    line_2 = gmsh.model.geo.addLine(p8, p4) # Curve 10

    # Create curve loops with the lines and arcs
    cl = gmsh.model.geo.addCurveLoop([-arc_i3, -arc_i2, -arc_i1, line_1, arc_e1, arc_e2, arc_e3, line_2])

    return cl
end

function generate_geometry(animal::RigidEllipse, h::Float64)
    # Define points based on animal attributes
    x1 = [animal.a*cos(animal.θ₀/2), animal.b*sin(animal.θ₀/2)]
    x2 = [0, animal.b]
    x3 = [-animal.a, 0]
    x4 = [0, -animal.b]
    x5 = [animal.a*cos(animal.θ₀/2), -animal.b*sin(animal.θ₀/2)]
    x6 = [(animal.a+animal.e)*cos(animal.θ₀/2), (animal.b+animal.e)*sin(animal.θ₀/2)]
    x7 = [0, animal.b+animal.e]
    x8 = [-(animal.a+animal.e), 0]
    x9 = [0, -(animal.b+animal.e)]
    x10 = [(animal.a+animal.e)*cos(animal.θ₀/2), -(animal.b+animal.e)*sin(animal.θ₀/2)]    

    # Transform points with rotation and translation
    aux1 = R(animal.α) * x1 .+ [animal.x, animal.y]
    aux2 = R(animal.α) * x2 .+ [animal.x, animal.y]
    aux3 = R(animal.α) * x3 .+ [animal.x, animal.y]
    aux4 = R(animal.α) * x4 .+ [animal.x, animal.y]
    aux5 = R(animal.α) * x5 .+ [animal.x, animal.y]
    aux6 = R(animal.α) * x6 .+ [animal.x, animal.y]
    aux7 = R(animal.α) * x7 .+ [animal.x, animal.y]
    aux8 = R(animal.α) * x8 .+ [animal.x, animal.y]
    aux9 = R(animal.α) * x9 .+ [animal.x, animal.y]
    aux10 = R(animal.α) * x10 .+ [animal.x, animal.y]
    

    # Add points in Gmsh
    cp = gmsh.model.geo.addPoint(animal.x, animal.y, 0, h)
    e_p1 = gmsh.model.geo.addPoint(aux1[1], aux1[2], 0, h)
    e_p2 = gmsh.model.geo.addPoint(aux2[1], aux2[2], 0, h)
    e_p3 = gmsh.model.geo.addPoint(aux3[1], aux3[2], 0, h)
    e_p4 = gmsh.model.geo.addPoint(aux4[1], aux4[2], 0, h)
    e_p5 = gmsh.model.geo.addPoint(aux5[1], aux5[2], 0, h)
    e_p6 = gmsh.model.geo.addPoint(aux6[1], aux6[2], 0, h)
    e_p7 = gmsh.model.geo.addPoint(aux7[1], aux7[2], 0, h)
    e_p8 = gmsh.model.geo.addPoint(aux8[1], aux8[2], 0, h)
    e_p9 = gmsh.model.geo.addPoint(aux9[1], aux9[2], 0, h)
    e_p10 = gmsh.model.geo.addPoint(aux10[1], aux10[2], 0, h)


    # Add lines and arcs in Gmsh for both the porous domain and the ellipse
    el1 = gmsh.model.geo.addEllipseArc(e_p1, cp, e_p3, e_p2)
    el2 = gmsh.model.geo.addEllipseArc(e_p2, cp, e_p3, e_p3)
    el3 = gmsh.model.geo.addEllipseArc(e_p3, cp, e_p3, e_p4)
    el4 = gmsh.model.geo.addEllipseArc(e_p4, cp, e_p3, e_p5)
    el5 = gmsh.model.geo.addEllipseArc(e_p6, cp, e_p8, e_p7)
    el6 = gmsh.model.geo.addEllipseArc(e_p7, cp, e_p8, e_p8)
    el7 = gmsh.model.geo.addEllipseArc(e_p8, cp, e_p8, e_p9)
    el8 = gmsh.model.geo.addEllipseArc(e_p9, cp, e_p8, e_p10)
    lel1 = gmsh.model.geo.addLine(e_p6, e_p1)
    lel2 = gmsh.model.geo.addLine(e_p5, e_p10)

    cl = gmsh.model.geo.addCurveLoop([el1, el2, el3, el4, lel2, -el8, -el7, -el6, -el5,  lel1])

    return cl, [el1, el2, el3, el4, el5, el6, el7, el8, lel1, lel2]
end

# Function to generate animals based on a definition
function generate_animals(definition)
    animals = Vector{MarineForms}()
    for (form_type, params) in definition
        for _ in 1:params[:N]
            flag = false
            while !flag
                new_animal = generate(form_type, params)
                if is_collision_free(new_animal, animals)
                    push!(animals, new_animal)
                    flag = true
                end
            end
        end
    end
    return animals
end

"This function checks if exists collision between two marine forms.
Since at this time we only have Clamps and Circles, we can use the exterior radius
to check the collision and the entry of the function can be an abstract animal. At the 
moment we introduce new animals, we can't use this directly, but we can use the concrete types."
# function check_collisions(form1::MarineForms, form2::MarineForms)
#     dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2)
    
#     return dist < (form1.re + form2.re + max(form1.tol, form2.tol) )  # Condition to have a collision
# end
function check_collisions(form1::Clamp, form2::Circle)
    dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2)
    
    return dist < (form1.re + form2.re + max(form1.tol, form2.tol) )  # Condition to have a collision
end
function check_collisions(form1::Clamp, form2::Clamp)
    dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2)
    
    return dist < (form1.re + form2.re + max(form1.tol, form2.tol) )  # Condition to have a collision
end

function check_collisions(form1::Circle, form2::Circle)
    dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2)
    
    return dist < (form1.re + form2.re + max(form1.tol, form2.tol) )  # Condition to have a collision
end
function check_collisions(form1::Clamp, form2::RigidEllipse)
    dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2)
    
    return dist < (form1.re + (form2.a + form2.e) + max(form1.tol, form2.tol) )  # Condition to have a collision
end

function check_collisions(form1::RigidEllipse, form2::Clamp)
    dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2)
    
    return dist < (form2.re + (form1.a + form1.e) + max(form2.tol, form1.tol) )  # Condition to have a collision
end

function check_collisions(form1::Circle, form2::RigidEllipse)
    dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2)
    
    return dist < (form1.re + (form2.a + form2.e) + max(form1.tol, form2.tol) )  # Condition to have a collision
end

function check_collisions(form1::RigidEllipse, form2::Circle)
    dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2)
    
    return dist < (form2.re + (form1.a + form1.e) + max(form1.tol, form2.tol) )  # Condition to have a collision
end

function check_collisions(form1::RigidEllipse, form2::RigidEllipse)
    dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2)
    
    return dist < ((form1.a + form1.e) + (form2.a + form2.e))  # Condition to have a collision
end

# + max(form1.tol, form2.tol)
"""
This function checks if exists collision between the animal created and the previously created animals.
Due to that the array initially is filled with nothing, the function checks if the animal is not nothing,
and then checks if the animal collides with the existing animals.
"""
function is_collision_free(new_animal::MarineForms, existing_animals::Vector{MarineForms})
    for animal in existing_animals
        if check_collisions(new_animal, animal)
            return false  # collision
        end
    end
    return true  
end

end 


