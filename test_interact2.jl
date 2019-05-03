using PlotlyJS
using Interact
using Blink

const strlen=20   #for the length of buttons

# any better way to save these things????
mutable struct Status
    clicked_points::Array{Tuple{Float32,Float32},1}
end
const status=Status(Array{Tuple{Float32,Float32},1}())



function genbuttons(contents,callbacks)
    buttons=Array{Widget{:button,Int64},1}()
    for (content,callback) in zip(contents,callbacks)
        if (length(content)>strlen) content=content[1:strlen] end
        nl=div(strlen-length(content),2)
        nr=strlen-nl-length(content)
        push!(buttons,button('_'^nl*content*'_'^nr))
        on(callback,buttons[end])
    end
    buttons
end


# plot the data in every column of yy
function maketraces(x,yy,labels=nothing)
    traces=Array{GenericTrace{Dict{Symbol,Any}},1}()
    for i=1:size(yy,2)
        if labels != nothing
            push!(traces,scatter(;x=x,y=yy[:,i],mode="lines",marker_color="black",name=labels[i]))
        else
            push!(traces,scatter(;x=x,y=yy[:,i],mode="lines",marker_color="black",marker_size=1))
        end
    end
    traces
end

# prepare the whole screen
function init_screen()
    layout=Layout(;hovermode="closest",showlegend=false,width=900,height=600)
    win=Window()
    return win,layout
end

global (win,layout)=init_screen();
global p1=plot(scatter(;x=0,y=0),layout);


# here's the part that actually makes the buttons
contents=Array{String,1}()
callbacks=Array{Any,1}()

# button for next
push!(contents,"next")
function next(nclick)
    traces=maketraces(rand(100),rand(100,10))
    react!(p1,traces,layout)
end
push!(callbacks,next)

# button for select window
push!(contents,"select window")
function selectwindow(ndum)
    relayout!(p1,Dict(:title=>"click on the trace to select window"))
    on(getcursor,p1["click"])
    on(addvertline,p1["click"])
end
function getcursor(data::Dict{String,Any})
    push!(status.clicked_points,
    (Float32(data["points"][1]["x"]),Float32(data["points"][1]["y"])))
end
function addvertline(data)
    if length(status.clicked_points)>2
        deletetraces!(p1,length(p1.plot.data)-1)
        popfirst!(status.clicked_points)
    end
    addtraces!(p1,scatter(;x=[status.clicked_points[end][1],status.clicked_points[end][1]],y=[0,1],marker_color="red"))
end
push!(callbacks,selectwindow)


# button for confirm and output xpos
push!(contents,"confirm window")
function confirmwindow(ndum)
    if (length(status.clicked_points) != 2)
        relayout!(p1,Dict(:title=>"something not right"))
    else
        relayout!(p1,Dict(:title=>"writing out $(status.clicked_points)"))
        deletetraces!(p1,length(p1.plot.data))
        deletetraces!(p1,length(p1.plot.data))
        pop!(status.clicked_points)
        pop!(status.clicked_points)
        off(p1["click"],getcursor)
        off(p1["click"],addvertline)
    end
end
push!(callbacks,confirmwindow)

# button for delete trace
push!(contents,"delete trace")
function deletetrace(ndum)
    relayout!(p1,Dict(:title=>"click on the trace to delete"))
    on(gettrace,p1["click"])
end
function gettrace(data)
    deletetraces!(p1,data["points"][1]["curveNumber"]+1)   #remember the plus 1
    off(p1["click"],gettrace)
end
push!(callbacks,deletetrace)

buttons=genbuttons(contents,callbacks)

body!(win,vbox(hbox(buttons),p1))
