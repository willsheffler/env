"""Tests for settings."""
from os import path

from EasyClangComplete.plugin.settings.settings_manager import SettingsManager
from EasyClangComplete.tests.gui_test_wrapper import GuiTestWrapper


class test_settings(GuiTestWrapper):
    """Test settings."""

    def test_setup_view(self):
        """Test that setup view correctly sets up the view."""
        file_name = path.join(path.dirname(__file__),
                              'test_files',
                              'test.cpp')
        self.check_view(file_name)
        self.tear_down()

    def test_init(self):
        """Test that settings are correctly initialized."""
        manager = SettingsManager()
        settings = manager.user_settings()
        self.assertIsNotNone(settings.verbose)
        self.assertIsNotNone(settings.include_file_folder)
        self.assertIsNotNone(settings.include_file_parent_folder)
        self.assertIsNotNone(settings.triggers)
        self.assertIsNotNone(settings.common_flags)
        self.assertIsNotNone(settings.clang_binary)
        self.assertIsNotNone(settings.flags_sources)
        self.assertIsNotNone(settings.errors_on_save)

    def test_valid(self):
        """Test validity."""
        manager = SettingsManager()
        settings = manager.user_settings()
        self.assertTrue(settings.is_valid())

    def test_populate_flags(self):
        """Testing include population."""
        # open any existing file
        file_name = path.join(path.dirname(__file__),
                              'test_files',
                              'test_wrong_triggers.cpp')
        self.set_up_view(file_name)
        # now test the things
        manager = SettingsManager()
        settings = manager.user_settings()
        self.assertTrue(settings.is_valid())

        initial_common_flags = list(settings.common_flags)
        settings = manager.settings_for_view(self.view)
        dirs = settings.common_flags

        current_folder = path.dirname(self.view.file_name())
        parent_folder = path.dirname(current_folder)
        self.assertTrue(len(initial_common_flags) <= len(dirs))
        self.assertTrue(initial_common_flags[0] in dirs)
        self.assertFalse(initial_common_flags[1] in dirs)
        self.assertTrue(("-I" + current_folder) in dirs)
        self.assertTrue(("-I" + parent_folder) in dirs)
        self.tear_down()
