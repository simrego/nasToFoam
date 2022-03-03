/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

Application
    nasToFoam

Group
    grpMeshConversionUtilities

Description
    nastran dat format mesh conversion.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "IFstream.H"
#include "StringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Dat file format.
enum struct FORMAT : char
{
    FREE = 0,       // Free format, column delimiter: ','
    SMALL = 8,      // Small format, every column 8 char wide
    LARGE = 16      // Large format, forst column 8, others 16 char wide.
};
// File format, small by default.
FORMAT format = FORMAT::SMALL;

// Last extracted entry kw
// Always read the next one if we are done with a line/multiline
string entryBuff;

// Buffer for last comment with it's line number. (patch, cellZone names)
word commentBuffer;
label commentLine;

// Get the next entry and store it in the buffer (fwd)
// Also return a reference to it.
string& getEntry(IFstream& is);

// Ignore everything until we find the "BEGIN BULK" entry
bool findBulk(IFstream& is)
{
    string line;
    while (is.good())
    {
        is.getLine(line);
        if (line.starts_with("BEGIN BULK")) return true;
    }
    return false;
}

// Process a commented line.
// e.g. NX nastran write the property card name as comment before the entry.
// So now just store the last word from the comment which is the name.
void processCommentedLine(IFstream& is)
{
    string line;
    is.getLine(line);
    // Yep, this is not too safe...
    size_t lastSpace = line.find_last_of(' ');
    if (lastSpace < line.size())
    {
        commentBuffer = line.substr(lastSpace + 1);
        commentLine = is.lineNumber();
    }
}

// Skip everything on this line, and next lines if we have a multiline entry.
// The only problem if we get the '\n' from the stream. In this case the
// next line will be skipped which shouldn't be...
// WARNING: Make sure the '\n' is kept in the stream. Only happens with 
// free format.
void finishEntry(IFstream& is)
{
    string line;
    while (true)
    {
        is.getLine(line);

        // Remove stupid carriage return char...
        line.erase(std::remove(line.begin(), line.end(), 13),
                 line.end());

        // If the end of the line is '+', the next line is the same entry.
        if (!line.ends_with('+')) break;
    }

}

// Extract the next column from the file.
string getColumn(IFstream& is, char size = char(format))
{
    if (size == 0) // Increase for free format to avoid overflow
    {
        size = 63;
    }
 
    char buff[63];  // More than enough.

    label i = 0;
    label nRead = 0;
    while (nRead++ < size && i < size)
    {
        is.get(buff[i]);
        // Stop at new line, or ',' for free format.
        if (buff[i] == '\n')
        {
            i++;
            break;
        }
        else if (buff[i] == ',') // only possible in free format
        {
            // Dont increase i, so we will simply change it to '\0'
            break;
        }
        else if (isspace(buff[i]) || buff[i] == 13)
        {
            // space or carriage return char (happens sometimes ?!?!?!)
            // Just continue without increasing i
            continue;
        }
        i++;
    }
    // Close the string.
    buff[i] = '\0';

    // Continue on new line if end with "+\n" and next start with + or *
    if (buff[0] == '+' && buff[i-1] == '\n' &&
        (is.peek() == '+' || is.peek() == '*'))
    {
        // Multiline, ignore the 1st column
        if (format == FORMAT::FREE)
        {
            char c;
            while (is.get(c) && c != ',') {}
        }
        else
        {
            is.readRaw(buff, 8);
        }

        // We are in the next line, start again.
        return getColumn(is, size);
    }

    // In free format '\n' is removed from the stream but we will need it
    // to properly finish the entry.
    // Put it back if we have it and close the string.
    if (format == FORMAT::FREE && buff[i-1] == '\n')
    {
        buff[i-1] = '\0';
        is.putback('\n');
    }

    return string(buff);
}

string& getEntry(IFstream& is)
{
    // Finish the current line/multiline
    finishEntry(is);
    // Process comments
    while (is.good() && is.peek() == '$') processCommentedLine(is);

    switch (format)
    {
    case FORMAT::SMALL:
    case FORMAT::LARGE:
        entryBuff = getColumn(is, 8);
        break;
    case FORMAT::FREE:
        entryBuff = getColumn(is, 0);
        break;
    }

    if (entryBuff.ends_with('*')) // We have multiline, ignore *.
        entryBuff.removeEnd('*');

    return entryBuff;
}

// Read the next column and return as label
label getLabel(IFstream& is)
{
    return readLabel(getColumn(is));
}

// Read the next column and return as scalar
// Scientific notation sucks in nastran...
// Sometimes we have E, sometimes we don't... ?!?!
scalar getScalar(IFstream& is)
{
    string data = getColumn(is);
    // This is stupid. And probably costly.
    label id = data.find('+', 1); // contains e+
    if (id == -1) id = data.find('-', 1); // contains e-
    if (id != -1 && (data[id-1] != 'e' && data[id-1] != 'E'))
    {
        data = data.substr(0, id) + "e" + data.substr(id);
    }

    return readFloat(data);
}

// Read points until we find some different entry.
// GRID card format, where CP is ignored:
// GRID   ID   CP   X  Y  Z  ...
void readPoints
(
    IFstream& is,
    DynamicList<point>& points,
    DynamicList<label>& pointIDs
)
{
    Info<< "\tReading \"GRID\" entries." << endl;

    label pointi;
    do
    {
        pointi = getLabel(is);
        if (pointi >= pointIDs.size())
        {
            pointIDs.resize(pointi + 1, -1);
        }
        pointIDs[pointi] = points.size();

        // Ignore CP column...
        getColumn(is);
        // Get the 3 coordinate
        point pt;
        pt[0] = getScalar(is);
        pt[1] = getScalar(is);
        pt[2] = getScalar(is);
        points.append(pt);

    } while (getEntry(is) == "GRID");
}

// Read cells of a given type until we find a different keyword.
template<cellModel::modelType TYPE>
void readCell
(
    string name,
    IFstream& is,
    DynamicList<cellShape>& cells,
    Map<DynamicList<label>>& cellPropIDs,
    DynamicList<label>& pointIDs
)
{
    Info<< "\tReading " << name << " entries." << endl;
    
    const cellModel& model = cellModel::ref(TYPE);
    label cellPropID;
    do
    {
        getColumn(is);   // ignore cell ID
        cellPropID = getLabel(is);
        if (!cellPropIDs.found(cellPropID))
        {
            cellPropIDs.insert(cellPropID, DynamicList<label>());
        }
        cellPropIDs.at(cellPropID).append(cells.size());
        
        labelList verts(model.nPoints());
        forAll(verts, i)
        {
            verts[i] = pointIDs[getLabel(is)];
        }
        cells.append(cellShape(model, verts, true));

    } while (getEntry(is) == name);
}

// Read faces until we find a different kieyword
template<label nPoints>
void readFaces
(
    string name,
    IFstream& is,
    Map<DynamicList<face>>& patches,
    DynamicList<label>& pointIDs
)
{
    label patchI;
    do
    {
        getColumn(is); // ignore ID
        patchI = getLabel(is);
        if (!patches.found(patchI))
        {
            patches.insert(patchI, DynamicList<face>());
        }

        face fVerts(nPoints);
        forAll(fVerts, i)
        {
            fVerts[i] = pointIDs[getLabel(is)];
        }
        patches.at(patchI).append(fVerts);

    } while (getEntry(is) == name);
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert nastran dat file to OpenFOAM. Units assumed to be in meters."
    );
    argList::noParallel();
    argList::addArgument(".dat file");
    argList::addOption(
        "format",
        "word",
        "Input format. Options: small(default), large, free"
    );
    argList::addBoolOption(
        "defaultNames",
        "Use default patch and cellZone names, don't use the comments."
    );

    #include "setRootCase.H"
    #include "createTime.H"

    if (args.found("format"))
    {
        word inpForm = args.getOrDefault<word>("format", "small");
        if (inpForm == "small") 
        {
            // Do nothing. Initialised to small.
        }
        else if (inpForm == "large")
        {
            format = FORMAT::LARGE;
        }
        else if (inpForm == "free")
        {
            format = FORMAT::FREE;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown format: " << inpForm << "."
                << exit(FatalError);
        }
    }
    bool defaultNames = args.found("defaultNames");

    const auto datName = args.get<fileName>(1);
    IFstream inFile(datName);

    if (!inFile.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << datName
            << exit(FatalError);
    }

    if (!findBulk(inFile))
    {
        FatalErrorInFunction
            << "Cannot find \"BEGIN BULK\" entry."
            << exit(FatalError);
    }

    // Points
    DynamicList<point> points;
    // Nastran indexing. pointIDs[nastranIndex] = <points index>
    // Could contain a lot of -1 elements, but the renumbering is
    // really fast in a cost of a little extra memory..
    DynamicList<label> pointIDs;
    // Cells
    DynamicList<cellShape> cells;
    // Cell property IDs
    // key is the nastran porperty ID
    Map<DynamicList<label>> cellPropIDs;
    // Patches
    // key is the nastran porperty ID
    Map<DynamicList<face>> patches;
    // Porperty card names
    Map<word> propNames;

    // Read the next entry into the buffer.
    getEntry(inFile);

    Info<< "Start reading file." << endl;
    while (inFile.good())
    {
        if (entryBuff == "GRID")
        {
            readPoints(inFile, points, pointIDs);
        }
        else if (entryBuff == "CTETRA")
        {
            readCell<cellModel::modelType::TET>(
                "CTETRA", inFile, cells, cellPropIDs, pointIDs
            );
        }
        else if (entryBuff == "CPYRAM")
        {
            readCell<cellModel::modelType::PYR>(
                "CPYRAM", inFile, cells, cellPropIDs, pointIDs
            );
        }
        else if (entryBuff == "CHEXA")
        {
            readCell<cellModel::modelType::HEX>(
                "CHEXA", inFile, cells, cellPropIDs, pointIDs
            );
        }
        else if (entryBuff == "CTRIA3")
        {
            readFaces<3>(
                "CTRIA3", inFile, patches, pointIDs
            );
        }
        else if (entryBuff == "CQUAD4")
        {
            readFaces<4>("CQUAD4", inFile, patches, pointIDs);
        }
        else if (entryBuff == "PSOLID" || entryBuff == "PSHELL")
        {
            // Property names
            label propI = getLabel(inFile);
            if (propNames.found(propI))
            {
                FatalErrorInFunction
                    << "Property ID: " << propI << " is already defined."
                    << exit(FatalError);
            }

            if (commentLine == inFile.lineNumber() && !defaultNames)
            {
                propNames.insert(propI, commentBuffer);
            }
            else
            {
                propNames.insert(propI, "");
            }
            getEntry(inFile);
        }
        else if (entryBuff == "ENDDATA")
        {
            Info<< "Finished reading file." << endl;
            break;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot process keyword: \"" << entryBuff.c_str()
                << "\", on line " << inFile.lineNumber() << "."
                << exit(FatalError);
        }
    }

    DynamicList<faceList> patchFaces;
    wordList patchNames;
    label unnamedPatchN = 0;
    forAll(patches, i)
    {
        label propI = patches.toc()[i];
        word propName = propNames.at(propI);
        if (patches[propI].size())
        {
            patchFaces.append(std::move(patches[propI]));
            if (propName.empty())
            {
                patchNames.append("patch_" + std::to_string(unnamedPatchN++));
            }
            else
            {
                patchNames.append(propName);
            }
        }
    }

    Info<< "Constructing the mesh." << endl;
    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        pointField(points),
        cells,
        patchFaces,
        patchNames,
        wordList(patchNames.size(), polyPatch::typeName),
        "defaultFaces",
        polyPatch::typeName,
        wordList()
    );

    
    if (cellPropIDs.size())
    {
        Info<< "Adding cell zones." << endl;

        List<cellZone*> cZones;
        label unnamedCellZoneN = 0;
        forAll(cellPropIDs, i)
        {
            label propI = cellPropIDs.toc()[i];
            word propName = propNames.at(propI);
            if (propName.empty())
                propName = "cellZone_" + std::to_string(unnamedCellZoneN++);
            cZones.append
            (
                new cellZone
                (
                    propName,
                    cellPropIDs[propI],
                    i,
                    mesh.cellZones()
                )
            );
        }

        mesh.addZones(List<pointZone*>(), List<faceZone*>(), cZones);
    }

    Info<< endl;
    Info<< "Mesh information:" << endl
        << "Number of points: " << mesh.nPoints() << endl
        << "Number of faces: " << mesh.nFaces() << endl
        << "Number of cells: " << mesh.nCells() << endl
        << "Patch names:" << endl;
    forAll(mesh.boundaryMesh(), i)
    {
        Info<< "\t" << mesh.boundaryMesh().get(i)->name() << endl;
    }
    if (mesh.cellZones().size())
    {
        Info<< "Cell zones:" << endl;
        forAll(mesh.cellZones(), i)
        {
            Info<< "\t" << mesh.cellZones()[i].name() << endl;
        }
    }
    Info<< endl;

    mesh.removeFiles();
    mesh.write();

    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
