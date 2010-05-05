#
# This file and all modifications and additions to the pristine
# package are under the same license as the package itself.
#

# norootforbuild

Name:           engrid
Version:	1.2rc1
Release:	1
Summary:	open-source mesh generator for CFD
Group:		Applications/Engineering
License:	GPL
Url:		http://engits.eu/engrid
Requires:	libqt4 vtk netgen
BuildRequires:  libqt4-devel vtk vtk-devel netgen netgen-devel
Source:		http://files.engits.eu/engrid-1.2rc1.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-build

%description

%prep
%setup

%build
cd src
export CXXFLAGS="$RPM_OPT_FLAGS"
export VTKINCDIR=%_includedir/vtk
export VTKLIBDIR=%_libdir
qmake
make

%install
mkdir -p %buildroot/%_bindir
%__cp engrid %buildroot/%_bindir

%clean
rm -fr %buildroot

%files
%defattr(-,root,root)
%doc engrid_manual.pdf
%_bindir/engrid

%changelog
