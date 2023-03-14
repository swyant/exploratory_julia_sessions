### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 3da59e96-c29a-11ed-237c-5f3e284ee15c
mutable struct Counter
	c::Int
end

# ╔═╡ 71e35dde-18d9-4677-99e3-42a4ff27ab13
function update!(c::Counter)
	c.c += 1
end

# ╔═╡ ecfbcdc4-12dc-4e0a-838f-5a6b93d0f510
cntr = Counter(0)

# ╔═╡ 4ebaecd9-c6cb-4bf3-a6ec-c2d7f137799b
begin 
	vedge = nothing
	update!(cntr)
end

# ╔═╡ 904cde90-3af5-4bcf-86a5-443e283709f0
cntr

# ╔═╡ a6e272bb-0cb3-46b4-aa95-e47c7ff48197
begin 
	vedge
	cntr
end

# ╔═╡ Cell order:
# ╠═3da59e96-c29a-11ed-237c-5f3e284ee15c
# ╠═71e35dde-18d9-4677-99e3-42a4ff27ab13
# ╠═ecfbcdc4-12dc-4e0a-838f-5a6b93d0f510
# ╠═4ebaecd9-c6cb-4bf3-a6ec-c2d7f137799b
# ╠═904cde90-3af5-4bcf-86a5-443e283709f0
# ╠═a6e272bb-0cb3-46b4-aa95-e47c7ff48197
