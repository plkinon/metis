import os

# Define the folder path
folder_path = "../classes/Integrator"

# Define the text to insert
text_to_insert = """            % Computes residual vector & tangent matrix
            %
            % :param zn1: state vector for next time step
            % :param zn: state vector at current time step
            % :param this_system: System object
            % :returns: [ResidualVector, TangentMatrix] for the Newton's method to update zn1

"""
text_lines = text_to_insert.splitlines()
Num_of_Lines_To_insert = len(text_lines)

# Flag the target function to be recognized
edit = "compute_resi_tang"


# Iterate through files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".m"):  # You can change the file extension to match your files
        file_path = os.path.join(folder_path, filename)

        # Read the lines
        with open(file_path, "r") as file:
            lines = file.readlines()

        # Search for the line that meets the specified condition
        for i, line in enumerate(lines):
            if (
                edit in line
                and "".join(lines[i + 1 : i + 1 + Num_of_Lines_To_insert])
                != text_to_insert
            ):
                delete = 1
                j = i + 1
                modified_lines = []
                while delete == 1:
                    if "%" in lines[j]:
                        j = j + 1
                    else:
                        delete = 0
                for n, line2 in enumerate(lines):
                    if n >= i + 1 and n < j + 1:
                        print("Removing line" + str(n))
                    else:
                        modified_lines.append(line2)

                lines = modified_lines

                print("Editing: " + filename)

                # Insert the new line after the matching line
                lines.insert(i + 1, text_to_insert)

        # Write the modified lines
        with open(file_path, "w") as file:
            file.writelines(lines)
