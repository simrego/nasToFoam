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
    FREE = 0,
    SMALL = 8,
    LARGE = 16
};
// File format, small by default.
FORMAT format = FORMAT::SMALL;

// Last extracted entry kw
string entryBuff;

// Buffer for last comment with it's line number. (patch, cellZone names)
word commentBuffer;
label commentLine;

// Get the next entry and store it in the buffer (fwd)
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
    label lastSpace = line.find_last_of(' ');
    if (lastSpace < line.size())
    {
        commentBuffer = line.substr(lastSpace + 1);
        commentLine = is.lineNumber();
    }
}

// Skip everything on this line, and next lines if we have a multiline entry.
// The only problem if we get the '\n' from the stream. In this case the
// next line will be skipped which shouldn't be...
// WARNING: Make sure the '\n' is kept in the stream.
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
 
    label i = 0;
    char buff[63];  // More than enough, get size no chars

    while (i < size)
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
            // DOnt increase i, so we will change it to '\0'
            break;
        }
        i++;
    }
    buff[i] = '\0';

    // Continue on new line if end with "+\n" and next start with + or *
    if (buff[0] == '+' && buff[i-1] == '\n' &&
        (is.peek() == '+' || is.peek() == '*'))
    {
        if (format == FORMAT::FREE)
        {
            char c;
            while (is.get(c) && c != ',') {}
        }
        else is.readRaw(buff, 8);    // Ignore 1st column

        // Do it again.
        return getColumn(is, size);
    }

    // In free format '\n' is removed from the stream but we will need it
    // to finish the entry. Put it back IF we have it.
    if (format == FORMAT::FREE && buff[i-1] == '\n')
    {
        buff[i-1] = '\0';
        is.putback('\n');
    }

    string data = string(buff);
    // Remove whitespaces
    data.erase(std::remove(data.begin(), data.end(), ' '),
                 data.end());

    // Remove stupid carriage return char...
    data.erase(std::remove(data.begin(), data.end(), 13),
             data.end());

    return data;
}

string& getEntry(IFstream& is)
{
    // Skip comments
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

label getLabel(IFstream& is)
{
    return readLabel(getColumn(is));
}

// Scientific notation sucks in nastran...
scalar getScalar(IFstream& is)
{
    string data = getColumn(is);
    // This is stupid. And probably costly.
    label id = data.find('+', 1); // contains e+
    if (id == -1) id = data.find('-', 1); // contains e-
    if (id != -1) data = data.substr(0, id) + "e" + data.substr(id);

    return readFloat(data);
}

// Read points until we find some different entry.
// GRID card format, where CP is ignored:
// GRID   ID   CP   X  Y  Z  ...
void readPoints(IFstream& is, DynamicList<point> &points)
{
    label pointi;
    while (true)
    {
        pointi = getLabel(is);
        if (pointi != points.size() + 1)
        {
            FatalErrorInFunction
                << "Points are not in order."
                << exit(FatalError);
        }

        // Ignore CP column...
        getColumn(is);
        // Get the 3 coordinate
        point pt;
        pt[0] = getScalar(is);
        pt[1] = getScalar(is);
        pt[2] = getScalar(is);
        points.append(pt);
        // Ignore the other entries
        finishEntry(is); //is.getLine(nullptr);

        // No more GRID entry.
        if (!getEntry(is).starts_with("GRID")) break;
    }
}

// Read cells of a given type.
template<cellModel::modelType TYPE>
void readCell
(
    string name,
    IFstream& is,
    DynamicList<cellShape>& cellVerts,
    DynamicList<label>& cellPropIDs
)
{
    const cellModel& model = cellModel::ref(TYPE);
    while (entryBuff == name)
    {
        // Ignore cell ID
        getLabel(is);
        cellPropIDs.append(getLabel(is));
        labelList verts(model.nPoints());
        forAll(verts, i)
        {
            // Vertex indexing starts from 1, but we store from 0.
            verts[i] = getLabel(is) - 1;
        }
        cellVerts.append(cellShape(model, verts, true));

        finishEntry(is);
        getEntry(is);
    }
}

template<label nPoints>
void readFaces
(
    string name,
    IFstream& is,
    DynamicList<face>& boundaryFaces,
    DynamicList<label>& facePropIDs
)
{
    while (entryBuff == name)
    {
        getLabel(is); // ignore ID
        facePropIDs.append(getLabel(is));
        face fVerts(nPoints);
        forAll(fVerts, i)
        {
            // Vertex indexing starts from 1, but we store from 0.
            fVerts[i] = getLabel(is) - 1;
        }
        boundaryFaces.append(fVerts);

        finishEntry(is); // Get the end of the line
        getEntry(is);
    }
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert nastran dat file to OpenFOAM. Only for small format and Meters."
    );
    argList::noParallel();
    argList::addArgument(".dat file");
    argList::addOption(
        "format",
        "word",
        "Input format. Options: small(default), large, free"
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
    // Cell vertices
    DynamicList<cellShape> cellVerts;
    // Cell property IDs
    DynamicList<label> cellPropIDs;
    // Boundary faces
    DynamicList<face> boundaryFaces;
    // Face property IDs
    DynamicList<label> facePropIDs;
    // Cell zoneIDs
    DynamicList<label> cellZoneIDs;
    wordList cellZoneNames;
    // Patch IDs
    DynamicList<label> patchIDs;
    wordList patchNames;

    // Read the next entry into the buffer.
    getEntry(inFile);

    Info << "Start reading file." << endl;
    while (inFile.good())
    {
        if (entryBuff == "GRID")
        {
            readPoints(inFile, points);
        }
        else if (entryBuff == "CTETRA")
        {
            readCell<cellModel::modelType::TET>(
                "CTETRA", inFile, cellVerts, cellPropIDs
            );
        }
        else if (entryBuff == "CPYRAM")
        {
            readCell<cellModel::modelType::PYR>(
                "CPYRAM", inFile, cellVerts, cellPropIDs
            );
        }
        else if (entryBuff == "CHEXA")
        {
            readCell<cellModel::modelType::HEX>(
                "CHEXA", inFile, cellVerts, cellPropIDs
            );
        }
        else if (entryBuff == "CTRIA3")
        {
            readFaces<3>("CTRIA3", inFile, boundaryFaces, facePropIDs);
        }
        else if (entryBuff == "PSOLID")
        {
            // Cell Zones
            if (commentLine == inFile.lineNumber())
            {
                cellZoneNames.append(commentBuffer);
            }
            else
            {
                cellZoneNames.append(
                    "cellZone" + std::to_string(cellZoneNames.size()));
            }
            cellZoneIDs.append(getLabel(inFile));
            finishEntry(inFile);
            getEntry(inFile);
        }
        else if (entryBuff == "PSHELL")
        {
            // Patches
            if (commentLine == inFile.lineNumber())
            {
                patchNames.append(commentBuffer);
            }
            else
            {
                patchNames.append(
                    "patch" + std::to_string(patchNames.size()));
            }
            patchIDs.append(getLabel(inFile));
            finishEntry(inFile);
            getEntry(inFile);
        }
        else if (entryBuff == "ENDDATA")
        {
            Info << "Finished reading file." << endl;
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
    patchFaces.resize(patchIDs.size());
    forAll(boundaryFaces, i)
    {
        label pI = patchIDs.find(facePropIDs[i]);
        patchFaces[pI].append(boundaryFaces[i]);
    }

    Info << "Constructing the mesh." << endl;
    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        pointField(points),
        cellVerts,
        patchFaces,
        patchNames,
        wordList(patchNames.size(), polyPatch::typeName),
        "defaultFaces",
        polyPatch::typeName,
        wordList()
    );

    
    if (cellZoneIDs.size())
    {
        Info << "Adding cell zones." << endl;

        List<cellZone*> cZones(cellZoneIDs.size());
        
        DynamicList<labelList> cellZones;
        cellZones.resize(cellZoneIDs.size());
        forAll(cellPropIDs, i)
        {
            label zI = cellZoneIDs.find(cellPropIDs[i]);
            cellZones[zI].append(i);
        }
        
        forAll(cellZoneIDs, i)
        {
            cZones[i] = new cellZone
            (
                cellZoneNames[i],
                cellZones[i],
                i,
                mesh.cellZones()
            );
        }

        mesh.addZones(List<pointZone*>(), List<faceZone*>(), cZones);
    }

    mesh.removeFiles();
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
