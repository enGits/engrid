#include "dialoglineedit.h"
#include "dialoglineeditplugin.h"

#include <QtPlugin>

DialogLineEditPlugin::DialogLineEditPlugin(QObject *parent)
    : QObject(parent)
{
  initialized = false;
}

void DialogLineEditPlugin::initialize(QDesignerFormEditorInterface * /* core */)
{
  if (initialized)
    return;

  initialized = true;
}

bool DialogLineEditPlugin::isInitialized() const
{
  return initialized;
}

QWidget *DialogLineEditPlugin::createWidget(QWidget *parent)
{
  return new DialogLineEdit(parent);
}

QString DialogLineEditPlugin::name() const
{
  return "DialogLineEdit";
}

QString DialogLineEditPlugin::group() const
{
  return "Tunnellicht";
}

QIcon DialogLineEditPlugin::icon() const
{
  return QIcon();
}

QString DialogLineEditPlugin::toolTip() const
{
  return "";
}

QString DialogLineEditPlugin::whatsThis() const
{
  return "";
}

bool DialogLineEditPlugin::isContainer() const
{
  return false;
}

QString DialogLineEditPlugin::domXml() const
{
  return "<widget class=\"DialogLineEdit\" name=\"dialogLineEdit\">\n"
         " <property name=\"geometry\">\n"
         "  <rect>\n"
         "   <x>0</x>\n"
         "   <y>0</y>\n"
         "   <width>100</width>\n"
         "   <height>25</height>\n"
         "  </rect>\n"
         " </property>\n"
         "</widget>\n";
}

QString DialogLineEditPlugin::includeFile() const
{
  return "dialoglineedit.h";
}

//! [0]
Q_EXPORT_PLUGIN2(dialoglineeditplugin, DialogLineEditPlugin)
//! [0]
