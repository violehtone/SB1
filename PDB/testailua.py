def assign_ss(phi, psi):
	""" Assign a secondary structure type based on the phi
	and psi angles of a residue """
	### START CODING HERE
	# for code checking purposes use the terms "loop", "alpha" or "beta"

	secondary_structure = ""

	if(psi > 0 and phi > 0):
		secondary_structure = "alpha" #left handed
	elif(psi > 0 and phi < 0):
		secondary_structure = "beta"
	elif(psi < 0 and phi > 0):
		secondary_structure = "loop"
	elif(psi < 0 and phi < 0):
		secondary_structure = "alpha" #right handed
    else:
        print("")

    return secondary_structure

print(assign_ss(60, 25))
